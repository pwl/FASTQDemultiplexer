"""

Extracts the sequences at positions `positions` from FASTQ.Records
`records` and saves them into a sequence `seq`

"""
function extract!(seq,records,positions)
    start = 1
    nrec = length(records)
    for i in 1:nrec
        if issubset(positions[i],1:length(records[i].sequence))
            range = records[i].sequence[positions[i]]
            copy!(seq,start,
                  records[i].data,first(range),length(range))
            start += length(range)
        else
            error("The sequence length of read $i is too short for the given protocol.")
        end
    end
    return seq
end


function recordview!(recout::FASTQ.Record,recin::FASTQ.Record,position::UnitRange)

    recout.data = recin.data
    recout.filled = recin.filled
    recout.identifier = recin.identifier
    recout.description = recin.description
    recout.sequence = recin.sequence[position]
    recout.quality = recin.quality[position]
    return recout
end


function recordview(record::FASTQ.Record,position::UnitRange)
    FASTQ.Record(record.data,
                 record.filled,
                 record.identifier,
                 record.description,
                 record.sequence[position],
                 record.quality[position])
end


function interpret!{N}(ir::InterpretedRecord{N},
                       interpreter::Interpreter{N})
    extract!(ir.cell, ir.records, interpreter.cellpos)
    extract!(ir.umi, ir.records, interpreter.umipos)
    # Bio.Seq.encode_copy!(ir.cellseq,1,ir.cell,1,length(ir.cell))

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    recordview!(ir.output, ir.records[insertid], insertpos)

    ir.cellid = gen_id(ir.cell)
    ir.umiid = gen_id(ir.umi)
    return nothing
end


function demultiplex{N}(ih::InputHandle{N},
                        interpreter::Interpreter{N},
                        outputs,
                        barcodes::Barcodes)

    ir = InterpretedRecord(interpreter)

    while !eof(ih)
        read!(ih, ir.records)

        # TODO: move the interpreter into InputHandles
        interpret!(ir, interpreter)

        if ir.umiid == 0 || ir.cellid == 0
            ir.unmatched = true
            ir.groupid, ir.groupname = -1, "unassigned"
        else
            ir.groupid, ir.groupname = getgroup(barcodes,ir.cellid)
            ir.unmatched = ir.groupid == -1
        end

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
type Demultiplexer{N}
    interpreter::Interpreter{N}
    barcodes::Barcodes
    # TODO: this is not type stable
    outputs::Vector
    inputdir::String
    # for debugging purposes
    maxreads
end


function demultiplex(dem::Demultiplexer)
    fastqfiles = listfastq(dem.inputdir,dem.interpreter)

    results = pmap(fastqfiles) do fastq
        input = InputHandle(fastq,maxreads=dem.maxreads)
        inputname = basename(input.name)
        println("Starting $inputname")

        outputs = map(dem.outputs) do out
            T, options = out
            tmpdir = joinpath(options[:outputdir],inputname)
            T(dem.interpreter; merge(options,Dict(:outputdir=>tmpdir))...)
        end

        demultiplex(input,dem.interpreter,outputs,dem.barcodes)
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
