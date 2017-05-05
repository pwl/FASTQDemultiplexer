__precompile__()
module FASTQDemultiplexer

using Bio
using Bio.Seq
using Bio.Seq.FASTQ

include("protocol.jl")
include("base.jl")
include("output.jl")
include("demultiplexer.jl")

end # module
