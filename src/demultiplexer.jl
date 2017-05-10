"""

Extracts the sequences at positions `positions` from FASTQ.Records
`records` and saves them into a sequence `seq`

"""
function extract!(seq,records,positions)
    start = 1
    nrec = length(records)
    for i in 1:nrec
        range = records[i].sequence[positions[i]]
        copy!(seq,start,
              records[i].data,first(range),length(range))
        start += length(range)
    end
    return seq
end


function interpret!{N}(ir::InterpretedRecord{N},
                       interpreter::Interpreter{N})
    extract!(ir.cell, ir.records, interpreter.cellpos)
    # Bio.Seq.encode_copy!(ir.cellseq,1,ir.cell,1,length(ir.cell))
    extract!(ir.umi, ir.records, interpreter.umipos)

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    # TODO: can we do without the copy here?
    ir.output = copy(ir.records[insertid])
    ir.output.sequence = ir.output.sequence[insertpos]
    ir.output.quality = ir.output.quality[insertpos]
    ir.cellid = gen_id(ir.cell)
    return nothing
end


function FASTQdemultiplex{N}(io::NTuple{N,IO},
                             interpreter::Interpreter{N},
                             oh::OutputHandler;
                             cellbarcodes=[],
                             kwargs...)
    readers = map(FASTQ.Reader,io)

    # TODO: move this out as a callback
    demultiplexcell = Demultiplexer(Vector{DNASequence}(cellbarcodes),
                                    n_max_errors=0)
    function acceptedbarcode(ir)
        if isempty(cellbarcodes)
            return true
        else
            return demultiplex(demultiplexcell,DNASequence(ir.cell))[1] > 0
        end
    end

    ir = InterpretedRecord(interpreter)

    while !any(map(eof,readers))
        for i in 1:N
            read!(readers[i],ir.records[i])
        end

        interpret!(ir,interpreter)
        if length(oh.selectedcells) > 0
            ir.unmatched = !(ir.cellid in oh.selectedcells)
        end

        write(oh,ir)
    end

    return nothing
end


function openfile(output,interpreter;ext=".fastq")
    unmatched = open(joinpath(output,"unmatched.fastq"),"w")
    function gendescryptor(cellid,accepted)
        if accepted
            filename = String(cellid)
            return open(joinpath(output,filename)*ext,"w")
        else
            return unmatched
        end
    end
end
