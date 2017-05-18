function demultiplex(yamlfile::String)

    config = YAML.load_file(yamlfile)

    interpreter = Interpreter(Symbol(get(config, "protocol", "none")))
    inputdir = get(config,"inputdir", ".")::String
    outputdir = get(config, "outputdir", joinpath(inputdir,"demultiplexed"))::String
    maxreads = get(config, "maxreads", Inf)
    towrite = map(Symbol,get(config, "output", ["insert","umi"])::Vector{String})
    barcodes = get(config, "cellbarcodes", "")::String
    jobs = get(config, "jobs", Sys.CPU_CORES)::Int
    mergeoutput = get(config, "merge", true)
    maxopenfiles = get(config, "maxopenfiles", 10)

    fastqfiles = listfastq(inputdir,interpreter)

    jobs = min(jobs,length(fastqfiles))-1
    if jobs > 0
        addprocs(jobs)
        @everywhere import FASTQDemultiplexer
    end

    tmpdir=joinpath(outputdir,"tmp")

    @sync @parallel for reads in fastqfiles
        hs = InputHandle(reads,maxreads=maxreads)
        println("Starting $(basename(hs.name))")
        output = OutputHandler(
            interpreter,
            outputdir = joinpath(tmpdir,hs.name),
            towrite = towrite,
            cellbarcodes = barcodes,
            maxopenfiles = maxopenfiles)
        FASTQdemultiplex(hs,interpreter,output)
        close(hs)
        close(output)
    end

    if mergeoutput
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
