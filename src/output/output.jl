abstract Output

function Output(sym::Symbol,interpreter::Interpreter; kwargs...)
    if sym == :split
        OutputSplit(interpreter; kwargs...)
    elseif sym == :filtered
        OutputFiltered(interpreter; kwargs...)
    else
        error("Unrecognized output type: $sym")
    end
end

include("filtered.jl")
include("split.jl")