using Bio
using Bio.Seq
using Bio.Seq.FASTQ
using Iterators
using BufferedStreams


supported_protocols()=["marsseq"]


immutable Interpreter{N}
    insertread::Int
    insertpos::UnitRange{Int}
    cellpos::NTuple{N,UnitRange{Int}}
    umipos::NTuple{N,UnitRange{Int}}
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
