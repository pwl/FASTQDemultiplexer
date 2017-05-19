####################
### Split output ###
####################


type OutputSplit{F}
    cellhandles::Dict{UInt,FASTQ.Writer}
    umihandles::Dict{UInt,IOStream}
    namegen::F
    maxopenfiles::Int
end


function OutputSplit(protocol::Interpreter;
                     outputdir::String = ".",
                     maxopenfiles::Int = 1000,
                     kwargs...)

    namegen(cellid, unmatched) =
        joinpath(outputdir,
                 unmatched ? "unmatched" : String(cellid))

    return OutputSplit(Dict{UInt,FASTQ.Writer}(),
                       Dict{UInt,IOStream}(),
                       namegen,
                       maxopenfiles)
end

"""

    Writes a read (`InterpretedRecord`) to corresponding output file(s).

    """
function Base.write(oh::OutputSplit, ir::InterpretedRecord)
    if !ir.unmatched
        closetoomanyfiles!(oh)
        write_insert(oh,ir)
        write_umi(oh,ir)
    end
end


"""

    Counts the number of open files and closes them if necessary.  Should
    be made more intelligent in the future.

    """
function closetoomanyfiles!{N}(oh::OutputSplit{N})
    if length(oh.cellhandles) + length(oh.umihandles) >= oh.maxopenfiles
        # close only the cell-like handles, leave the raw handles open
        map(close, values(oh.cellhandles))
        map(close, values(oh.umihandles))

        oh.cellhandles=Dict()
        oh.umihandles=Dict()
    end
    return nothing
end


function write_insert(oh::OutputSplit,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid

    celldesc = get!(oh.cellhandles, writeid) do
        filename = oh.namegen(ir.cell, ir.unmatched)*".fastq"
        FASTQ.Writer(open(filename, "w+"))
    end
    write(celldesc,ir.output)

    return nothing
end


function write_umi(oh::OutputSplit,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid
    umidesc = get!(oh.umihandles, writeid) do
        open(oh.namegen(ir.cell, ir.unmatched)*".umi", "w+")
    end
    unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
    write(umidesc,'\n')

    return nothing
end


function Base.close(oh::OutputSplit)
    map(close, values(oh.cellhandles))
    map(close, values(oh.umihandles))
end
