abstract type Output end

# mergeoutput{O<:Output}(::Vector{O}; kwargs...) = nothing
Base.close(::Output) = nothing

include("filtered.jl")
include("split.jl")
