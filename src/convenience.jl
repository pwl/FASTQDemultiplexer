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

    inputs = map(fastqfiles) do reads
        InputHandle(reads,maxreads=config[:maxreads])
    end

    symtocons = Dict("split" => OutputSplit,
                     "filtered"=> OutputFiltered)

    cons = map(config[:output]) do out
        @show out
        symtocons[out[:type]], out
    end

    demultiplexer = Demultiplexer(inputs,interpreter,cellbarcodes,cons)

    FASTQdemultiplex(demultiplexer)

    return nothing

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
