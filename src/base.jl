type InterpretedRecord{N}
    umi::Vector{UInt8}
    cell::Vector{UInt8}
    cellseq::DNASequence
    output::FASTQ.Record
    records::NTuple{N,FASTQ.Record}
    unmatched::Bool
    cellid::UInt
    umiid::UInt
end


function InterpretedRecord{N}(interpreter::Interpreter{N})
    umi = Array(UInt8,sum(map(length,interpreter.umipos)))
    cell = Array(UInt8,sum(map(length,interpreter.cellpos)))
    cellseq = DNASequence(sum(map(length,interpreter.cellpos)))
    output = FASTQ.Record()
    records = map(x->FASTQ.Record(),1:N)
    return InterpretedRecord{N}(umi,cell,cellseq,output,(records...),false,0,0)
end
