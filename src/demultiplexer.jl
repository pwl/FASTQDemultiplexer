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
    # Bio.Seq.encode_copy!(ir.cellseq,1,ir.cell,1,length(ir.cell))
    if writeumis
        extract!(ir.umi, recs, interpreter.umipos)
    end

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    ir.output = recs[insertid]
    ir.output.sequence = ir.output.sequence[insertpos]
    ir.output.quality = ir.output.quality[insertpos]
    return nothing
end


function FASTQdemultiplex{N}(io::NTuple{N,IO},
                             interpreter::Interpreter{N};
                             output::String="output",
                             closebuffers::Bool=true,
                             writeumis::Bool=false,
                             cellbarcodes=[],
                             celldescryptor=openfile(output,interpreter,ext=".fastq"),
                             umidescryptor=openfile(output,interpreter,ext=".umi"),
                             kwargs...)
    readers = map(FASTQ.Reader,io)
    records = map(x->FASTQ.Record(),io)

    cellhandles = Dict{Vector{UInt8},FASTQ.Writer}()
    # TODO: IO is an abstract type, can this intruduce type instability?
    #umihandles = Dict{Int,IO}()

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

    function accepted(cellid)

    end
    ir = InterpretedRecord(interpreter)

    while !any(map(eof,readers))
        for i in 1:N
            read!(readers[i],records[i])
        end

        interpret!(ir,records,interpreter,writeumis)

        # TODO: have a custom write function
        celldesc = get!(cellhandles,ir.cell) do
            FASTQ.Writer(celldescryptor(ir.cell,acceptedbarcode(ir)))
        end
        write(celldesc,ir.output)

        # if writeumis
        #     umidesc = get!(umihandles,cellid) do
        #         umidescryptor(cellid)
        #     end
        #     unsafe_write(umidesc,pointer(ir.umi),length(ir.umi))
        #     write(umidesc,'\n')
        # end
    end

    if closebuffers
        map(close,values(cellhandles))
        # map(close,values(umihandles))
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
