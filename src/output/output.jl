abstract Output

# mergeoutput{O<:Output}(::Vector{O}; kwargs...) = nothing
Base.close(::Output) = nothing

include("filtered.jl")
include("split.jl")
