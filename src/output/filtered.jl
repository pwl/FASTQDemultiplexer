#######################
### Filtered output ###
#######################


type OutputFiltered{N} <: Output
    handles::NTuple{N,FASTQ.Writer}
end


function OutputFiltered{N}(protocol::Interpreter{N};
                           outputdir::String = ".",
                           kwargs...)

    handles = map(protocol.readnames) do name
        filename = joinpath(outputdir, name*".fastq")
        FASTQ.Writer(open(filename, "w"))
    end

    OutputFiltered{N}(handles)
end


Base.close(o::OutputFiltered) = map(close,o.handles)


function Base.write{N}(oh::OutputFiltered{N},ir::InterpretedRecord{N})
    if ! ir.unmatched
        for i in 1:N
            write(oh.handles[i],ir.records[i])
        end
    end
end
