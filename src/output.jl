type OutputHandler{F}
    cellhandles::Dict{Vector{UInt8},Vector{FASTQ.Writer}}
    umihandles::Dict{Vector{UInt8},IO}
    namegen::F
    writeumi::Bool
    # writeraw::Bool
    writeunmatched::Bool
    # outputdir::String
end


function OutputHandler(yamlfile)
    config = YAML.load_file(yamlfile)

    namegen(cellid)=String(cellid)

    writeumi = haskey(config,"writeumi") ? config["writeumi"] : false
    writeunmatched = haskey(config,"writeunmatched") ? config["writeunmatched"] : true

    return OutputHandler(Dict(), Dict(), namegen, writeumi, writeunmatched)
end


function Base.write(oh::OutputHandler, ir::InterpretedRecord)
    if !ir.unmatched

        celldesc = get!(oh.cellhandles, ir.cell) do
            filename = oh.namegen(ir.cell)*".fastq"
            FASTQ.Writer(open(filename, "w"))
        end
        write(celldesc,ir.output)

        if oh.writeumis

            umidesc = get!(oh.umihandles,cellid) do
                oh.namegen(cellid)*".umi"
            end
            unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
            write(umidesc,'\n')

        end

    end
    return nothing
end
