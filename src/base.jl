type InterpretedRecord
    umi::Vector{UInt8}
    cell::Vector{UInt8}
    cellseq::DNASequence
    output::FASTQ.Record
    unmatched::Bool
    cellid::UInt
end


function InterpretedRecord(interpreter::Interpreter)
    umi = Array(UInt8,sum(map(length,interpreter.umipos)))
    cell = Array(UInt8,sum(map(length,interpreter.cellpos)))
    cellseq = DNASequence(sum(map(length,interpreter.cellpos)))
    output = FASTQ.Record()
    return InterpretedRecord(umi,cell,cellseq,output,false,0)
end
