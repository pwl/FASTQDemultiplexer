using FASTQDemultiplexer
using Bio.Seq.FASTQ
using Base.Test

import FASTQDemultiplexer: idtoname, gen_id, Barcodes, Barcode,
    Interpreter, InterpretedRecord, getgroup, interpret!

# utils
@testset "utils" begin
    ain = UInt8['A','C','T','A','C','A','T','G','A','C']
    id = 0x00000000000184b2
    @test gen_id(ain) == id
    @test idtoname(id,length(ain))==String(ain)
end

# barcodes
@testset "barcodes" begin
    ain1 = UInt8['A','C','T','A','C','A','T','G','A','C']
    ain2 = reverse(ain1)
    ain3 = UInt8['C','A','A','A','A','A','T','G','A','C']
    ain4 = reverse(ain3)
    ainnone = UInt8['A','A','A','A','A','A','A','A','A','A']

    group1 = Set(map(gen_id,[ain1, ain2]))
    group2 = Set(map(gen_id,[ain3, ain4]))

    barcodes = Barcodes((group1,group2),("group1","group2"))

    @test typeof(Barcode(ain1)) == Barcode{length(ain1)}

    @test getgroup(barcodes,Barcode(ain1)) == (1, "group1")
    @test getgroup(barcodes,Barcode(ain4)) == (2, "group2")
    @test getgroup(barcodes,Barcode(ainnone)) == (-1, "unassigned")

    barcodes2 = Barcodes(["data/barcodes.txt"],["some-barcodes"])
    @test getgroup(barcodes2,Barcode(Vector{UInt8}("AAGG"))) == (1,"some-barcodes")
end


@testset "interpreter" begin

    mars = Interpreter(:marsseq)
    drop = Interpreter(:dropseq)
    x10 = Interpreter(:x10)

    @test typeof(mars) == Interpreter{2,10,6}
    @test typeof(drop) == Interpreter{2,12,8}
    @test typeof(x10) == Interpreter{3,16,10}

    @test mars.insertread == 1
    @test mars.insertpos == 8:60
    @test mars.cellpos == (4:7,1:6)
    @test mars.umipos == (1:0,7:12)

    @test drop.insertread == 2
    @test drop.insertpos == 1:39
    @test drop.cellpos == (1:12,1:0)
    @test drop.umipos == (13:20,1:0)

    @test x10.insertread  == 2
    @test x10.insertpos == 1:98
    @test x10.cellpos == (1:16,1:0,1:0)
    @test x10.umipos == (17:26,1:0,1:0)

end

@testset "record" begin
    N = 2; C = 4; U = 2;
    inter = Interpreter{N,C,U}(1,5:12,(2:3,3:4),(1:0,1:2),("R1","R2"))
    ir = InterpretedRecord(inter)

    ir.records = map(("data/r1.fastq","data/r2.fastq")) do file
        open(FASTQ.Reader,file) do f
            read(f)
        end
    end

    @test typeof(ir) == InterpretedRecord{N,C,U}
    @test length(ir.umi) == U
    @test length(ir.cell) == C

    cell = zeros(UInt8,C)
    FASTQDemultiplexer.extract!(cell,ir.records,inter.cellpos)
    @test cell == Vector{UInt8}("AAGG")

    umi = zeros(UInt8,U)
    FASTQDemultiplexer.extract!(umi,ir.records,inter.umipos)
    @test umi == Vector{UInt8}("TT")

    bcodes = Barcodes(["data/barcodes.txt"],["some-barcodes"])

    interpret!(ir,inter,bcodes)

    @test ir.cell == cell
    @test ir.umi == umi
    @test ir.cellid == Barcode(cell)
    @test ir.umiid == Barcode(umi)
    @test ir.groupid == 1
    @test ir.groupname == "some-barcodes"


    # read the second record
    ir.records = map(("data/r1.fastq","data/r2.fastq")) do file
        open(FASTQ.Reader,file) do f
            read(f)
            read(f)
        end
    end

    cell = zeros(UInt8,C)
    FASTQDemultiplexer.extract!(cell,ir.records,inter.cellpos)
    @test cell == Vector{UInt8}("CCTT")

    umi = zeros(UInt8,U)
    FASTQDemultiplexer.extract!(umi,ir.records,inter.umipos)
    @test umi == Vector{UInt8}("GG")

    interpret!(ir,inter,bcodes)

    @test ir.cell == cell
    @test ir.umi == umi
    @test ir.cellid == Barcode(cell)
    @test ir.umiid == Barcode(umi)
    @test ir.groupid == -1
    @test ir.groupname == "unassigned"
end
