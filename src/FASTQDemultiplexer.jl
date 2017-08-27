module FASTQDemultiplexer

using BioSequences.FASTQ
using YAML

include("utils.jl")
include("protocol.jl")
include("barcodes.jl")
include("record.jl")
include("input.jl")
include("output/output.jl")
include("demultiplexer.jl")

end # module
