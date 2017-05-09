using YAML
using Iterators

immutable Interpreter{N}
    insertread::Int
    insertpos::UnitRange{Int}
    cellpos::NTuple{N,UnitRange{Int}}
    umipos::NTuple{N,UnitRange{Int}}
    readnames::NTuple{N,String}
end


function Interpreter(yamlfile::String)
    config = YAML.load_file(yamlfile)

    nreads = length(config["reads"])

    cellpos = Array(UnitRange{Int}, nreads)
    umipos  = Array(UnitRange{Int}, nreads)
    readnames = Array(String, nreads)
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

    end

    return Interpreter{nreads}(insertread,insertpos,
                               (cellpos...),
                               (umipos...),
                               (readnames...))

end
