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

    outputs = pmap(fastqfiles) do reads
        hs = InputHandle(reads,maxreads=config[:maxreads])
        println("Starting $(basename(hs.name))")
        output = generateoutput(interpreter,config,hs.name)
        FASTQdemultiplex(hs,interpreter,output,cellbarcodes)
        close(hs)
        map(close,output)
        output
    end

    for i in 1:length(config[:output])
        println("Merging $(config[:output][i][:type])")
        outputi = [out[i] for out in outputs]
        @time mergeoutput(outputi; outputdir = config[:output][i][:outputdir])
    end

    return nothing

end


function generateoutput(interpreter,config,subdir)
    output = map(config[:output]) do out

        outputdir = joinpath(out[:outputdir],subdir)

        # TODO: move to Output?
        Output(Symbol(out[:type]),interpreter;
               merge(out,Dict(:outputdir=>outputdir))...)
        # TODO: should this work too?
        # Output(Symbol(out[:type]),interpreter;
        #        out..., outputdir = outputdir)
    end
    # TODO: is a tuple faster then a vector?
    return (output...)
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
