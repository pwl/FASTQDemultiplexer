####################
### Split output ###
####################


type OutputSplit{F} <: Output
    cellhandles::Dict{UInt,FASTQ.Writer}
    umihandles::Dict{UInt,IOStream}
    cellpaths::Dict{UInt,String}
    namegen::F
    maxopenfiles::Int
    outputdir::String
end


function OutputSplit(protocol::Protocol;
                     outputdir::String = ".",
                     maxopenfiles::Int = 1000,
                     kwargs...)

    mkpath(outputdir)

    namegen(cellid, unmatched) =
        joinpath(outputdir,
                 unmatched ? "unmatched" : String(cellid))

    return OutputSplit(Dict{UInt,FASTQ.Writer}(),
                       Dict{UInt,IOStream}(),
                       Dict{UInt,String}(),
                       namegen,
                       maxopenfiles,
                       outputdir)
end

"""

    Writes a read (`InterpretedRecord`) to corresponding output file(s).

    """
function Base.write(oh::OutputSplit, ir::InterpretedRecord)
    if !ir.unmatched
        closetoomanyfiles!(oh)
        writeinsert(oh,ir)
        writeumi(oh,ir)
        addpath(oh,ir)
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


function writeinsert(oh::OutputSplit,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid.val

    celldesc = get!(oh.cellhandles, writeid) do
        filename = oh.namegen(ir.cell, ir.unmatched)*".fastq"
        FASTQ.Writer(open(filename, "w+"))
    end
    write(celldesc,ir.output)

    return nothing
end


function writeumi(oh::OutputSplit,ir::InterpretedRecord)

    writeid = ir.unmatched ? UInt(0) : ir.cellid.val
    umidesc = get!(oh.umihandles, writeid) do
        open(oh.namegen(ir.cell, ir.unmatched)*".umi", "w+")
    end
    unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
    write(umidesc,'\n')

    return nothing
end


function addpath(oh::OutputSplit,ir::InterpretedRecord)
    writeid = ir.unmatched ? UInt(0) : ir.cellid.val
    get!(oh.cellpaths, writeid) do
        oh.namegen(ir.cell, ir.unmatched)
    end
    return nothing
end


function mergeoutput{N}(outputs::Vector{OutputSplit{N}};
                        outputdir::String = ".",
                        kwargs...)

    # 1) merge the cell barcodes
    paths = Dict{UInt,Vector{String}}()
    for o in outputs
        for (id,path) in o.cellpaths
            pbase = get!(paths, id) do
                String[]
            end
            push!(pbase,path)
        end
    end

    # 2) for each cellid concatenate the contents of all files into
    # one and remove the directories
    pmap(values(paths)) do ps
        mergecellid(outputdir,ps,".umi")
        mergecellid(outputdir,ps,".fastq")
    end

    # 3) clean up the directories
    pmap(outputs) do oh
        rm(oh.outputdir,force=true,recursive=true)
    end
end


function mergecellid(outputdir,path,ext)
    outfile = joinpath(outputdir,basename(path[1])*ext)
    infiles = map(f->f*ext,path)
    catfiles(outfile,infiles)
    map(rm, infiles)
end


function Base.close(oh::OutputSplit)
    map(close, values(oh.cellhandles))
    map(close, values(oh.umihandles))
end
