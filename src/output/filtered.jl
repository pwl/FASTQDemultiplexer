#######################
### Filtered output ###
#######################


type OutputFiltered{N} <: Output
    handles::NTuple{N,FASTQ.Writer}
    filenames::NTuple{N,String}
    outputdir::String
end


function OutputFiltered{N}(protocol::Interpreter{N};
                           outputdir::String = ".",
                           kwargs...)

    mkpath(outputdir)

    filenames = map(protocol.readnames) do name
        joinpath(outputdir, name*".fastq")
    end

    handles = map(filenames) do filename
        FASTQ.Writer(open(filename, "w+"))
    end

    OutputFiltered{N}(handles,filenames,outputdir)
end


Base.close(o::OutputFiltered) = map(close,o.handles)


function Base.write{N}(oh::OutputFiltered{N},ir::InterpretedRecord{N})
    if ! ir.unmatched
        for i in 1:N
            write(oh.handles[i],ir.records[i])
        end
    end
end


function mergeoutput{N}(outputs::Vector{OutputFiltered{N}};
                        outputdir::String = ".",
                        kwargs...)
    pmap(1:N) do i
        name = joinpath(outputdir,"$i.fastq")
        filenames = [oh.filenames[i] for oh in outputs]
        catfiles(name,filenames)
        map(rm,filenames)
    end

    pmap(outputs) do oh
        rm(oh.outputdir,force=true,recursive=true)
    end
end
