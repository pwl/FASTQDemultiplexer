type InterpretedRecord
    umi::Vector{UInt8}
    cell::Vector{UInt8}
    cellseq::DNASequence
    output::FASTQ.Record
end


function InterpretedRecord(interpreter::Interpreter)
    umi = Array(UInt8,sum(map(length,interpreter.umipos)))
    cell = Array(UInt8,sum(map(length,interpreter.cellpos)))
    cellseq = DNASequence(sum(map(length,interpreter.cellpos)))
    output = FASTQ.Record()
    return InterpretedRecord(umi,cell,cellseq,output)
end

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


function interpret!{N}(ir::InterpretedRecord,
                       recs::NTuple{N,FASTQ.Record},
                       interpreter::Interpreter{N},
                       writeumis)
    extract!(ir.cell, recs, interpreter.cellpos)
    Bio.Seq.encode_copy!(ir.cellseq,1,ir.cell,1,length(ir.cell))
    if writeumis
        extract!(ir.umi, recs, interpreter.umipos)
    end

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    # TODO how slow is the copy here?
    ir.output = recs[insertid]
    ir.output.sequence = ir.output.sequence[insertpos]
    ir.output.quality = ir.output.quality[insertpos]
    return nothing
end


function FASTQdemultiplex{N}(io::NTuple{N,IO},interpreter::Interpreter{N};
                             output="output",
                             closebuffers=true,
                             writeumis=false,
                             celldescryptor=openfile(output,interpreter,ext=".fastq"),
                             umidescryptor=openfile(output,interpreter,ext=".umi"),
                             kwargs...)
    readers = map(FASTQ.Reader,io)
    records = map(x->FASTQ.Record(),io)

    cellhandles = Dict{Int,FASTQ.Writer}()
    # TODO: IO is an abstract type, can this intruduce type instability?
    umihandles = Dict{Int,IO}()
    demultiplexcell = Demultiplexer(interpreter.cellbarcodes,n_max_errors=0)
    demultiplexumi  = Demultiplexer(interpreter.umibarcodes, n_max_errors=0)
    ir = InterpretedRecord(interpreter)

    while !any(map(eof,readers))
        for i in 1:N
            read!(readers[i],records[i])
        end

        interpret!(ir,records,interpreter,writeumis)
        cellid, _ = demultiplex(demultiplexcell,ir.cellseq)

        # TODO: have a custom write function
        celldesc = get!(cellhandles,cellid) do
            FASTQ.Writer(celldescryptor(cellid))
        end
        write(celldesc,ir.output)

        if writeumis
            umidesc = get!(umihandles,cellid) do
                umidescryptor(cellid)
            end
            unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
            write(umidesc,'\n')
        end
    end

    if closebuffers
        map(close,values(cellhandles))
        map(close,values(umihandles))
    end

    return nothing
end


function openfile(output,interpreter;ext=".fastq")
    unmatched = open(joinpath(output,"unmatched.fastq"),"w")
    function gendescryptor(cellid)
        if cellid > 0
            filename = string(interpreter.cellbarcodes[cellid])
            return open(joinpath(output,filename)*ext,"w")
        else
            return unmatched
        end
    end
end
