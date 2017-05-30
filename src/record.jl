type InterpretedRecord{N,C,U}
    umi::Vector{UInt8}
    cell::Vector{UInt8}
    output::FASTQ.Record
    records::NTuple{N,FASTQ.Record}
    unmatched::Bool
    cellid::Barcode{C}
    umiid::Barcode{U}
    groupid::Int
    groupname::String
end


function InterpretedRecord{N,C,U}(interpreter::Interpreter{N,C,U})
    umi = Array(UInt8,U)
    cell = Array(UInt8,C)
    output = FASTQ.Record()
    records = map(x->FASTQ.Record(),1:N)
    return InterpretedRecord{N,C,U}(
        umi,cell,output,(records...),false,Barcode{C}(0),Barcode{U}(0),-1,"")
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


"""

Extracts the sequences at positions `positions` from FASTQ.Records
`records` and saves them into a sequence `seq`

"""
function extract!(seq::Vector{UInt8},records,positions)
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


function interpret!{N,C,U}(ir::InterpretedRecord{N,C,U},
                           interpreter::Interpreter{N,C,U},
                           barcodes::Barcodes)

    extract!(ir.cell, ir.records, interpreter.cellpos)
    extract!(ir.umi, ir.records, interpreter.umipos)

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    recordview!(ir.output, ir.records[insertid], insertpos)

    ir.cellid = Barcode{C}(gen_id(ir.cell))
    ir.umiid = Barcode{U}(gen_id(ir.umi))

    if !isnull(ir.umiid) || !isnull(ir.cellid)
        ir.groupid, ir.groupname = -1, "unassigned"
        ir.unmatched = true
    else
        ir.groupid, ir.groupname = getgroup(barcodes,ir.cellid)
        ir.unmatched = ir.groupid == -1
    end

    return nothing
end
