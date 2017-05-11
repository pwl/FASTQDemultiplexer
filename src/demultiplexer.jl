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


function FASTQdemultiplex{N}(ih::InputHandle{N},
                             interpreter::Interpreter{N},
                             oh::OutputHandler)

    ir = InterpretedRecord(interpreter)

    while !eof(ih)
        read!(ih, ir.records)

        # TODO: move the interpreter into InputHandles
        interpret!(ir, interpreter)

        if ir.umiid == 0 || ir.cellid == 0
            ir.unmatched = true
        elseif length(oh.selectedcells) > 0
            ir.unmatched = !(ir.cellid in oh.selectedcells)
        end

        write(oh,ir)
    end

    return nothing
end


function demultiplex(yamlfile::String)

    config = YAML.load_file(yamlfile)

    interpreter = Interpreter(Symbol(get(config, "protocol", "none")))
    inputdir = get(config,"inputdir", ".")::String
    outputdir = get(config, "outputdir", joinpath(inputdir,"demultiplexed"))::String
    maxreads = get(config, "maxreads", Inf)
    towrite = map(Symbol,get(config, "output", ["insert","umi"])::Vector{String})
    barcodes = get(config, "cellbarcodes", "")::String
    jobs = get(config, "jobs", Sys.CPU_CORES)::Int

    input = FASTQDemultiplexer.InputHandler(inputdir,interpreter,maxreads=maxreads)
    reads = collect(take(input.handles,1))
    # reads = input.handles

    rm(outputdir,force=true,recursive=true)
    mkdir(outputdir)

    if jobs > 0
        addprocs(min(jobs,length(reads)))
        @everywhere import FASTQDemultiplexer
    end

    @time @sync @parallel for hs in reads
        println("Starting $(basename(hs.name))")
        output = OutputHandler(interpreter,
                               outputdir = joinpath(outputdir,hs.name),
                               towrite = towrite,
                               cellbarcodes = barcodes)
        FASTQdemultiplex(hs,interpreter,output)
        close(output)
    end
    return nothing
end
