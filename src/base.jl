using Bio
using Bio.Seq
using Bio.Seq.FASTQ
using Iterators


supported_protocols()=["marsseq"]


immutable Interpreter{N}
    insertread::Int
    insertpos::Range{Int}
    cellpos::NTuple{N,Range{Int}}
    umipos::NTuple{N,Range{Int}}
    cellbarcodes::Vector{DNASequence}
    umibarcodes::Vector{DNASequence}
end

"""
Generate all the combinations of ACGT of length `len`
"""
generatebarcodes(len) = collect(map(v->DNASequence(*(v...)),product(repeated(["A","C","G","T"],len)...)))


function Interpreter(id)
    if id == :marsseq
        subcellbc = DNASequence["CTACCA","CATGCT","GCACAT","TGCTCG","AGCAAT","AGTTGC","CCAGTT","TTGAGC","ACCAAC","GGTCCA","GTATAA","TTCGCT","AACTTG","CACATC","TCGGAA","AAGGAT","GCACAC","TCTGGC","CATAGC","CAGGAG","TGTCGG","ATTATG","CCTACC","TACTTA","GAAGAA","AGGATC","GACAGT","CCTATG","TCGCCT","ATAGCG","GTCGCC","ATTCTA","CGTTAC","GTCTGA","TTACGC","TTGAAT","AAGACA","CAGCAA","TCCAGC","CCAGAG","TCCTTG","AGGTTA","GTCATC","CCTTCG","TCTCGG","ATTGTC","GAACCT","TAATGA","AACGCA","CAACTC","CTGTAA","TAAGCA","AATGTT","CAGCGG","GACCAG","TATCCA","ACAGGT","CCAACA","GCCGTC","TATCTG","ACTAAG","CGCCTT","GCCTAG","TGCTGC","AGGTAA","CTAACT","GTAACA","TGTAAT","ATTATC","CTATGC","GTCCAC","TTGTCT","AACAAT","ATTCCT","GAAGGA","TCGCTA","ACAGTT","CAATAG","GACCGT","TCTGCA","ACTGTA","CCAGCA","GACCTA","TGCAAG","AGCATG","CGCTAT","GATATC","TGTAAC","AGGTCG","CTGCGG","GCAGCC","TTAATC","AGGTGC","CTGTGG","GCCGCA","TTATAT","CGGAAT","TGCAGC","CATTGA","CTGATG","AACTGG","TCCAGT","ATCGTC","TAGAGC","AAGGCT","GTAGCA","GCGATA","TAGTCG","GGTACC","GTATCG","GTACTA","GAACTG","GACTGA","AGCATC","TACGTA","CTAGTG","ATTGGC","CAGGTT","GCTATG","ATCTCG","ATACGC","CGACGT","TTCCGA","GTCTCA","ATCCGT","AGCTGT","CATCGT","CGCAGT","TGAGAC","ACGATG","TAGACT","TCACAG","TAACCG","GATCAC","ACGTAC","GGACTT","TCGTAG","CTGTAC","TTAGCC","GAGCTC","CTTGAC","GTAGTC","TGTACA","CCAAGT","CATCAG","GTCCAA","CAGTCA","TGCTGA","AGCTTA","CGTAGA","TCACGT","AGACTC","TAAGCT","TCGCGA","GCGTCA","AATCGG","TTCAGG","CTATAG","TCAAGC","GAGTCT","CCGTAA","TAGGCA","CGAGAT","TCAATG","CGAATG","ACTGGA","GGTTCA","TCAGTC","AAGCTT","GCTCTA","TAGCTA","TTAACG","TCAGCA","ACTCTG","CGTACG","GTGCAC","AGTACT","AGTCAG","CCTAGG","ACTAGT","AGTCTA","GCGTAT","CTCAGA","AGCGCT","GTCAAG","TAGCGT","ACGGTC","GCATTG","GGCTAA","CTGTGA","CATGCA","GATCGA"]
        poolbc = DNASequence["AGTC","CATG","TTGG","ACAG","CTAC","ATCA","TGAT","TCTA"]
        cellbc = collect(map(join,product(poolbc,subcellbc)))
        umibc = generatebarcodes(6)

        #     1  4   8           60
        # R1: XXXCCCCIIIIII---IIII
        #     1     7
        # R2: CCCCCCUUUUUU
        return Interpreter{2}(
            1,8:60,
            (4:7,1:6),
            (1:0,7:12),
            cellbc,
            umibc
        )
    end
end


type InterpretedRecord{S<:Bio.Seq.BioSequence}
    umi::S
    cell::S
    output::FASTQ.Record
end


function InterpretedRecord(interpreter::Interpreter)
    umi = DNASequence(sum(map(length,interpreter.umipos)))
    cell = DNASequence(sum(map(length,interpreter.cellpos)))
    output = FASTQ.Record()
    return InterpretedRecord(umi,cell,output)
end

"""

Extracts the sequences at positions `positions` from FASTQ.Records
`records` and saves them into a sequence `seq`

"""
function extract!(seq,records,positions)
    start = 1
    for (record,position) in zip(records,positions)
        copy!(seq,start,
              Bio.Seq.FASTQ.sequence(record,position),1)
        start += length(position)
    end
    return seq
end


function interpret!{N}(ir::InterpretedRecord,
                       recs::NTuple{N,FASTQ.Record},
                       interpreter::Interpreter{N};
                       outputquality = true,
                       writeumis = false,
                       kwargs...)
    extract!(ir.cell, recs, interpreter.cellpos)
    if writeumis
        extract!(ir.umi, recs, interpreter.umipos)
    end

    insertid = interpreter.insertread
    insertpos = interpreter.insertpos

    # TODO how slow is the copy here?
    ir.output = copy(recs[insertid])
    ir.output.sequence = recs[insertid].sequence[insertpos]
    if outputquality
        ir.output.quality = recs[insertid].quality[insertpos]
    else
        ir.output.quality = 1:0
    end
    return ir
end


function FASTQdemultiplex{N}(io::NTuple{N,IO},interpreter;
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
        for (reader,record) in zip(readers,records)
            read!(reader,record)
        end

        interpret!(ir,records,interpreter; writeumis = writeumis, kwargs...)
        cellid, _ = demultiplex(demultiplexcell,ir.cell)

        # TODO: have a custom write function
        celldesc = get!(cellhandles,cellid) do
            FASTQ.Writer(celldescryptor(cellid))
        end
        write(celldesc,ir.output)

        if writeumis
            umidesc = get!(umihandles,cellid) do
                umidescryptor(cellid)
            end
            write(umidesc,string(ir.umi)*"\n")
            # for (rec,rng) in  zip(records,interpreter.umipos)
            #     unsafe_write(umidesc,pointer(rec.data)+first(rec.sequence[rng])-1,length(rng))
            # end
            # write(umidesc,"\n")
        end
    end

    if closebuffers
        map(close,values(cellhandles))
        map(close,values(umihandles))
    end

    return nothing
end


function openfile(output,interpreter;ext=".fastq")
    function gendescryptor(cellid)
        if cellid > 0
            filename = string(interpreter.cellbarcodes[cellid])
        else
            filename = "unmatched"
        end
        return open(joinpath(output,filename)*ext,"w")
    end
end
