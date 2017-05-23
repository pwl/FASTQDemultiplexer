using Weave


type OutputQC <: Output
    cellid::Vector{UInt}
    umiid::Vector{UInt}
    cellbclen::UInt
    umibclen::UInt
    selectedcells::Vector{UInt}
end


function OutputQC(protocol::Interpreter;
                  outputdir::String = ".",
                  selectedcells::String = "",
                  kwargs...)
    bcodes = genselectedcells(selectedcells)
    cellbclen = sum(length(pos) for pos in protocol.cellpos)
    umibclen  = sum(length(pos) for pos in protocol.umipos)

    return OutputQC([],[],cellbclen,umibclen,bcodes)
end


function Base.write(o::OutputQC,ir::InterpretedRecord)
    push!(o.cellid,ir.cellid)
    push!(o.umiid,ir.umiid)
end


function mergeoutput(outputs::Vector{OutputQC};
                     outputdir::String = ".",
                     kwargs...)

    o1 = outputs[1]
    args = Dict{String,Any}()
    args["counts"], args["cellids"], args["umiids"] = counttable(outputs)
    args["cellnames"] = map(i->idtoname(i,o1.cellbclen),args["cellids"])
    args["uminames"]  = map(i->idtoname(i,o1.umibclen), args["umiids"])

    weave(Pkg.dir("FASTQDemultiplexer","src","output","QC","counts.jmd"),
          args = args, out_path = outputdir)
end


function counttable(outputs::Vector{OutputQC})
    #1) merge the encountered cells and umis
    cellids = Set{UInt}()
    umiids = Set{UInt}()
    @time for o in outputs
        union!(cellids,o.cellid)
        union!(umiids,o.umiid)
    end
    delete!(cellids,0)
    delete!(umiids,0)
    @time cellind = Dict(id=>i for (i,id) in enumerate(cellids))
    @time umiind  = Dict(id=>i for (i,id) in enumerate(umiids))

    @time counts = zeros(UInt,(length(cellids),length(umiids)))
    @time for o in outputs
        for (cellid, umiid) in zip(o.cellid, o.umiid)
            if cellid > 0 && umiid > 0
                counts[cellind[cellid],umiind[umiid]] += 1
            end
        end
    end
    counts, collect(cellids), collect(umiids)
end
