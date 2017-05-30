function demultiplex{N}(ih::InputHandle{N},
                        protocol::Protocol{N},
                        outputs,
                        barcodes::Barcodes)

    ir = InterpretedRecord(protocol)

    while !eof(ih)
        read!(ih, ir.records)

        # TODO: move the protocol into InputHandles
        interpret!(ir, protocol, barcodes)

        for oh in outputs
            write(oh,ir)
        end
    end

    return nothing
end


# TODO this type is a duck tape work and needs to be improved.  To do
# that I have to improve the input type for it to use an open--like
# constructor: store the file names in one type and open them to
# construct another (closable) type.
type Demultiplexer{N,C,U}
    protocol::Protocol{N,C,U}
    barcodes::Barcodes
    # TODO: this is not type stable
    outputs::Vector
    inputdir::String
    # for debugging purposes
    maxreads
end


function demultiplex(dem::Demultiplexer)
    fastqfiles = listfastq(dem.inputdir,dem.protocol)

    results = pmap(fastqfiles) do fastq
        input = InputHandle(fastq,maxreads=dem.maxreads)
        inputname = basename(input.name)
        println("Starting $inputname")

        outputs = map(dem.outputs) do out
            T, options = out
            tmpdir = joinpath(options[:outputdir],inputname)
            T(dem.protocol; merge(options,Dict(:outputdir=>tmpdir))...)
        end

        demultiplex(input,dem.protocol,outputs,dem.barcodes)
        close(input)
        map(close,outputs)
        outputs
    end

    for (i,(res, out)) in enumerate(zip(results,dem.outputs))
        println("Merging $(eltype(res))")
        resi = [r[i] for r in results]
        @time mergeoutput(resi; outputdir = out[2][:outputdir])
    end

    return nothing
end
