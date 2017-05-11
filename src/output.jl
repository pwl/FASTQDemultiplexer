"""

A type designed to handle the writing of the demultiplexed reads to
corresponding FASTQ files.

"""

type OutputHandler{N,F}
    cellhandles::Dict{UInt,FASTQ.Writer}
    umihandles::Dict{UInt,IOStream}
    rawhandles::NTuple{N,FASTQ.Writer}
    namegen::F
    towrite::Vector{Symbol}
    maxopenfiles::Int
    # todo move from here?
    selectedcells::Vector{UInt}
end


"""

Generates the OutputHandler given a path to a YAML file.

"""
function OutputHandler{N}(protocol::Interpreter{N};
                          cellbarcodes::String = "",
                          outputdir::String = ".",
                          maxopenfiles::Int = 1000,
                          towrite::Vector{Symbol} = [])

    supported = [:umi,:raw,:insert,:quality,:unmatched]
    unsupported = setdiff(towrite,supported)
    if ! isempty(unsupported)
        error("Unsupported symbols: $unsupported, use one of $supported")
    end

    selectedcells=UInt[]
    if cellbarcodes != ""
        if isfile(cellbarcodes)
            selectedcells = map(readdlm(cellbarcodes,String)[:,1]) do b
                @assert length(b) == sum(map(length,protocol.cellpos))
                gen_id(Vector{UInt8}(b))
            end
        else
            error("Could not find the file $bc")
        end
    end


    if !isdirpath(outputdir)
        mkdir(outputdir)
    end


    outputinsert = joinpath(outputdir,"insert")
    if !isdirpath(outputinsert)
        mkdir(outputinsert)
    end


    function namegen(cellid, unmatched)
        if unmatched
            basename = "unmatched"
        else
            basename = String(cellid)
        end
        return joinpath(outputinsert,basename)
    end


    if :raw in towrite
        outputraw = joinpath(outputdir,"raw")
        if !isdirpath(outputraw)
            mkdir(outputraw)
        end
        rawhandles = map(protocol.readnames) do name
            filename = joinpath(outputraw,name*".fastq")
            FASTQ.Writer(open(filename, "w+"), quality_header = :quality in towrite)
        end
    else
        rawhandles = map(protocol.readnames) do name
            FASTQ.Writer(DevNull)
        end
    end


    return OutputHandler(Dict{UInt,FASTQ.Writer}(),
                         Dict{UInt,IOStream}(),
                         (rawhandles...),
                         namegen,
                         towrite,
                         maxopenfiles, selectedcells)
end

"""

Writes a read (`InterpretedRecord`) to corresponding output file(s).

"""
function Base.write(oh::OutputHandler, ir::InterpretedRecord)

    if ir.unmatched && ! (:unmatched in oh.towrite)
        return nothing
    else
        closetoomanyfiles!(oh)

        if :raw in oh.towrite
            write_raw(oh,ir)
        end

        if :insert in oh.towrite
            write_insert(oh,ir)
        end

        if :umi in oh.towrite
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


"""

Writes the insert to the respective sink.

"""
function write_insert(oh::OutputHandler,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid

    celldesc = get!(oh.cellhandles, writeid) do
        filename = oh.namegen(ir.cell, ir.unmatched)*".fastq"
        FASTQ.Writer(open(filename, "w+"), quality_header = :quality in oh.towrite)
    end
    write(celldesc,ir.output)

    return nothing
end


"""

Writes a UMI to the UMI sink.

"""
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
