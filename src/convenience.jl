const defvalues = Dict(:inputdir => ".",
                       :outputdir => "output",
                       :maxreads => Inf,
                       :output => [Dict(:type => "split")],
                       :jobs => Sys.CPU_CORES,
                       :mergeoutput => false,
                       :protocol => "none",
                       :cellbarcodes => "")


keystosym(d)=d
keystosym(d::AbstractArray)=map(keystosym,d)

function keystosym(d::Dict)
    Dict(Symbol(k)=>keystosym(v) for (k,v) in d)
end

function demultiplex(yamlfile::String)

    userconfig = keystosym(YAML.load_file(yamlfile))
    config = merge(defvalues,userconfig)

    interpreter = Interpreter(Symbol(config[:protocol]))

    cellbarcodes = genselectedcells(config[:cellbarcodes])
    fastqfiles = listfastq(config[:inputdir],interpreter)

    jobs = min(config[:jobs],length(fastqfiles))-1
    if jobs > 0
        addprocs(jobs)
        @everywhere import FASTQDemultiplexer
    end

    @sync @parallel for reads in fastqfiles
        hs = InputHandle(reads,maxreads=config[:maxreads])
        println("Starting $(basename(hs.name))")
        output = generateoutput(interpreter,config,hs.name)
        FASTQdemultiplex(hs,interpreter,output,cellbarcodes)
        close(hs)
        map(close,output)
    end

    if config[:merge]
        println("Merging files...")
        mergeall(tmpdir,outputdir)
        rm(tmpdir,force=true,recursive=true)
    end

    return nothing

end


function mergeall(outputdir,saveto)
    flist=String[]
    rawlist=String[]
    for dir in readdir(outputdir)
        dirraw = joinpath(outputdir,dir,"raw")
        if isdir(dirraw)
            for f in readdir(dirraw)
                push!(rawlist,joinpath(dirraw,f))
            end
        end

        dirinsert = joinpath(outputdir,dir,"insert")
        if isdir(dirinsert)
            for f in readdir(dirinsert)
                push!(flist,joinpath(dirinsert,f))
            end
        end
    end

    mergebarcodes(flist,"umi",joinpath(saveto,"insert"))
    mergebarcodes(flist,"fastq",joinpath(saveto,"insert"))
    mergebarcodes(rawlist,"fastq",joinpath(saveto,"raw"))
end


function mergebarcodes(flist,ext,saveto)
    mkpath(saveto)

    flistfiltered = filter(f->ismatch(Regex("."*ext*"\$"),basename(f)),flist)
    basenames = map(basename,flistfiltered)

    @sync @parallel for bn in collect(take(unique(basenames),500))
        positions = basenames .== bn
        files = flistfiltered[positions]
        outputfile = joinpath(saveto,bn)
        catfiles(outputfile,files)
    end
    return nothing
end


function catfiles(output,files)
    fout = open(output,"w")
    for f in files
        write(fout,read(f))
    end
    close(fout)
end


function genselectedcells(cellbarcodes::String)
    if cellbarcodes != ""
        if isfile(cellbarcodes)
            selectedcells = map(readdlm(cellbarcodes,String)[:,1]) do b
                gen_id(Vector{UInt8}(b))
            end
            return Set{UInt}(selectedcells)
        else
            error("Could not find the file $bc")
        end
    else
        return selectedcells=Set{UInt}[]
    end
end


function generateoutput(interpreter,config,subdir)
    output = map(config[:output]) do out

        outputdir = joinpath(out[:outputdir],subdir)

        # TODO: move to Output?
        mkpath(outputdir)
        Output(Symbol(out[:type]),interpreter;
               merge(out,Dict(:outputdir=>outputdir))...)
        # TODO: should this work too?
        # Output(Symbol(out[:type]),interpreter;
        #        out..., outputdir = outputdir)
    end
    # TODO: is a tuple faster then a vector?
    return (output...)
end
