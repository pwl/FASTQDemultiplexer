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
    # Bio.Seq.encode_copy!(ir.cellseq,1,ir.cell,1,length(ir.cell))
    extract!(ir.umi, ir.records, interpreter.umipos)

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    recordview!(ir.output, ir.records[insertid], insertpos)

    ir.cellid = gen_id(ir.cell)
    return nothing
end


function Base.read!{N}(readers::NTuple{N,FASTQ.Reader},ir::InterpretedRecord{N})
    for i in 1:N
        read!(readers[i],ir.records[i])
    end
end


function FASTQdemultiplex{N}(readers::NTuple{N,FASTQ.Reader},
                             interpreter::Interpreter{N},
                             oh::OutputHandler;
                             kwargs...)

    ir = InterpretedRecord(interpreter)

    while !any(map(eof,readers))
        read!(readers,ir)

        # TODO: move the interpreter into InputHandles
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
