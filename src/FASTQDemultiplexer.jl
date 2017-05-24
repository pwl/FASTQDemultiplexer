module FASTQDemultiplexer

using Bio
using Bio.Seq
using Bio.Seq.FASTQ

include("utils.jl")
include("protocol.jl")
include("input.jl")
include("base.jl")
include("barcodes.jl")
include("output/output.jl")
include("demultiplexer.jl")

end # module
