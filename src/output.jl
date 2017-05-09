type OutputHandler{F}
    cellhandles::Dict{UInt,FASTQ.Writer}
    umihandles::Dict{UInt,IOStream}
    namegen::F
    writeumi::Bool
    # writeraw::Bool
    writeunmatched::Bool
    maxopenfiles::Int
    selectedcells::Vector{UInt}
end


function OutputHandler(yamlfile)
    config = YAML.load_file(yamlfile)

    writeumi = haskey(config,"writeumi") ? config["writeumi"] : false
    writeunmatched = haskey(config,"writeunmatched") ? config["writeunmatched"] : false
    output = haskey(config,"output") ? config["output"] : "."
    maxopenfiles = haskey(config,"maxopenfiles") ? config["maxopenfiles"] : 100
    cellbarcodes= haskey(config,"cellbarcodes") ? config["cellbarcodes"] : String[]
    selectedcells = map(gen_id, map(DNASequence,cellbarcodes))

    function namegen(cellid, unmatched)
        if unmatched
            basename = "unmatched"
        else
            basename = String(cellid)
        end
        return joinpath(output,basename)
    end

    return OutputHandler(Dict{UInt,FASTQ.Writer}(),
                         Dict{UInt,IOStream}(),
                         namegen,
                         writeumi, writeunmatched,
                         maxopenfiles, selectedcells)
end


function Base.write(oh::OutputHandler, ir::InterpretedRecord)

    if ir.unmatched && ! oh.writeunmatched
        return nothing
    end

    if length(oh.cellhandles) + length(oh.umihandles) >= oh.maxopenfiles
        error("Too many open files")
        close(oh)
        oh.cellhandles=Dict()
        oh.umihandles=Dict()
    end

    writeid = ir.unmatched ? UInt(0) : ir.cellid

    celldesc = get!(oh.cellhandles, writeid) do
        filename = oh.namegen(ir.cell, ir.unmatched)*".fastq"
        FASTQ.Writer(open(filename, "w+"))
    end
    write(celldesc,ir.output)

    if oh.writeumi
        umidesc = get!(oh.umihandles, writeid) do
            open(oh.namegen(ir.cell, ir.unmatched)*".umi", "w+")
        end
        unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
        write(umidesc,'\n')
    end

    return nothing
end

function Base.close(oh::OutputHandler)
    map(close, values(oh.cellhandles))
    map(close, values(oh.umihandles))
end
