immutable Protocol{N,C,U}
    insertread::Int
    insertpos::UnitRange{Int}
    cellpos::NTuple{N,UnitRange{Int}}
    umipos::NTuple{N,UnitRange{Int}}
    readnames::NTuple{N,String}
    readlengths::NTuple{N,Int}
end

# const config_path = joinpath(Pkg.dir("FASTQDemultiplexer"),"config","protocols")
const config_path = joinpath(@__DIR__,"..","config","protocols")
const protocol_config_paths = Dict(:marsseq => joinpath(config_path,"marsseq.yml"),
                                   :dropseq => joinpath(config_path,"dropseq.yml"),
                                   :x10 => joinpath(config_path,"10x.yml"))

function Protocol(prot::Symbol)

    allowed_protocols = keys(protocol_config_paths)

    if prot in allowed_protocols
        return Protocol(protocol_config_paths[prot])
    else
        error("Unknown protocol: $prot, use one of $allowed_protocol")
    end
end


function Protocol(yamlfile::String)
    config = YAML.load_file(yamlfile)

    nreads = length(config["reads"])

    cellpos = Array{UnitRange{Int}}(nreads)
    umipos  = Array{UnitRange{Int}}(nreads)
    readnames = Array{String}(nreads)
    readlengths = Array{Int}(nreads)
    insertpos = 1:0
    insertread = 0
    insertdefined = false

    # parse the reads
    for (i,read) in enumerate(config["reads"])

        if haskey(read,"cell")
            # TODO add safer parsing?
            cellpos[i] = eval(parse(read["cell"]))
        else
            cellpos[i] = 1:0
        end

        if haskey(read,"umi")
            umipos[i] = eval(parse(read["umi"]))
        else
            umipos[i] = 1:0
        end

        if haskey(read,"insert")
            if insertdefined
                error("The field insert can be defined in only one read")
            end
            insertpos = eval(parse(read["insert"]))
            insertread = i
            insertdefined = true
        end

        if haskey(read,"name")
            readnames[i]=read["name"]
        else
            error("Name of the read is required")
        end

        if haskey(read,"length")
            readlengths[i]=read["length"]
        else
            error("Length of the read is required")
        end

    end

    C = sum(map(length,cellpos))
    U = sum(map(length,umipos))

    return Protocol{nreads,C,U}(insertread,
                                insertpos,
                                (cellpos...),
                                (umipos...),
                                (readnames...),
                                (readlengths...))

end

cellidtoname(id::UInt,::Protocol{N,C,U}) where {N,C,U} = idtoname(id,C)

umiidtoname(id::UInt,::Protocol{N,C,U}) where {N,C,U} = idtoname(id,U)
