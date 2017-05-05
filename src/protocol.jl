using YAML
using Iterators

supported_protocols()=["marsseq"]


immutable Interpreter{N}
    insertread::Int
    insertpos::UnitRange{Int}
    cellpos::NTuple{N,UnitRange{Int}}
    umipos::NTuple{N,UnitRange{Int}}
    cellbarcodes::Vector{DNASequence}
    umibarcodes::Vector{DNASequence}
end


function Interpreter(yamlfile::String)
    config = YAML.load_file(yamlfile)
    try
        nreads = length(config["reads"])

        cellpos = Array(UnitRange{Int}, nreads)
        umipos  = Array(UnitRange{Int}, nreads)
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

        end

        cellbarcodes = DNASequence[]
        umibarcodes = DNASequence[]

        if haskey(config,"barcodes") && config["barcodes"] != nothing
            barcodes = config["barcodes"]

            # parse the barcodes
            if haskey(barcodes,"cells")
                cellbarcodes = DNASequence[bc for bc in barcodes["cells"]]
            end

            if haskey(barcodes,"umis")
                umibarcodes = DNASequence[bc for bc in barcodes["umis"]]
            end
        end

        return Interpreter{nreads}(insertread,insertpos,
                                   (cellpos...),
                                   (umipos...),
                                   cellbarcodes,
                                   umibarcodes)

    catch err
        warn("Unable to parse the YAML file.")
        throw(err)
    end
end
