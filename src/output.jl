"""

A type designed to handle the writing of the demultiplexed reads to
corresponding FASTQ files.

"""

type OutputHandler{N,F}
    cellhandles::Dict{UInt,FASTQ.Writer}
    umihandles::Dict{UInt,IOStream}
    rawhandles::NTuple{N,FASTQ.Writer}
    namegen::F
    writeumi::Bool
    writeraw::Bool
    writetrimmed::Bool
    writeunmatched::Bool
    maxopenfiles::Int
    # todo move from here?
    selectedcells::Vector{UInt}
end

"""

Generates the OutputHandler given a path to a YAML file.

"""
function OutputHandler{N}(protocol::Interpreter{N}, yamlfile; subdirectory = "")
    config = YAML.load_file(yamlfile)

    writeumi = haskey(config,"writeumi") ? config["writeumi"] : false
    writeraw = haskey(config,"writeraw") ? config["writeraw"] : false
    writetrimmed = haskey(config,"writetrimmed") ? config["writetrimmed"] : false
    writeunmatched = haskey(config,"writeunmatched") ? config["writeunmatched"] : false

    maxopenfiles = haskey(config,"maxopenfiles") ? config["maxopenfiles"] : 100
    cellbarcodes= haskey(config,"cellbarcodes") ? config["cellbarcodes"] : String[]
    selectedcells = map(gen_id, map(DNASequence,cellbarcodes))

    output = haskey(config,"output") ? config["output"] : "."
    output = joinpath(output,subdirectory)
    if !isdirpath(output)
        mkdir(output)
    end

    outputtrimmed = joinpath(output,"trimmed")
    if !isdirpath(outputtrimmed)
        mkdir(outputtrimmed)
    end

    function namegen(cellid, unmatched)
        if unmatched
            basename = "unmatched"
        else
            basename = String(cellid)
        end
        return joinpath(outputtrimmed,basename)
    end


    if writeraw
        outputraw = joinpath(output,"raw")
        if !isdirpath(outputraw)
            mkdir(outputraw)
        end
        rawhandles = map(protocol.readnames) do name
            filename = joinpath(outputraw,name*".fastq")
            FASTQ.Writer(open(filename, "w+"))
        end
    end

    return OutputHandler(Dict{UInt,FASTQ.Writer}(),
                         Dict{UInt,IOStream}(),
                         (rawhandles...),
                         namegen,
                         writeumi, writeraw, writetrimmed, writeunmatched,
                         maxopenfiles, selectedcells)
end

"""

Writes a read (`InterpretedRecord`) to corresponding output file(s).

"""
function Base.write(oh::OutputHandler, ir::InterpretedRecord)

    if ir.unmatched && ! oh.writeunmatched
        return nothing
    else
        closetoomanyfiles!(oh)

        if oh.writeraw
            write_raw(oh,ir)
        end

        if oh.writetrimmed
            write_trimmed(oh,ir)
        end

        if oh.writeumi
            write_umi(oh,ir)
        end
    end
end


"""

Counts the number of open files and closes them if necessary.  Should
be made more intelligent in the future.

"""
function closetoomanyfiles!{N}(oh::OutputHandler{N})
    if length(oh.cellhandles) + length(oh.umihandles) >= oh.maxopenfiles
        # close only the cell-like handles, leave the raw handles open
        map(close, values(oh.cellhandles))
        map(close, values(oh.umihandles))

        oh.cellhandles=Dict()
        oh.umihandles=Dict()
    end
    return nothing
end


function write_raw{N}(oh::OutputHandler{N},ir::InterpretedRecord{N})
    writeid = ir.unmatched ? UInt(0) : ir.cellid

    for i in 1:N
        write(oh.rawhandles[i],ir.records[i])
    end

end


function write_trimmed(oh::OutputHandler,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid

    celldesc = get!(oh.cellhandles, writeid) do
        filename = oh.namegen(ir.cell, ir.unmatched)*".fastq"
        FASTQ.Writer(open(filename, "w+"))
    end
    write(celldesc,ir.output)

    return nothing
end


function write_umi(oh::OutputHandler,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid
    umidesc = get!(oh.umihandles, writeid) do
        open(oh.namegen(ir.cell, ir.unmatched)*".umi", "w+")
    end
    unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
    write(umidesc,'\n')

    return nothing
end


function Base.close(oh::OutputHandler)
    map(close, values(oh.cellhandles))
    map(close, values(oh.umihandles))
    map(close, oh.rawhandles)
end
