using FASTQDemultiplexer
using Bio.Seq
using Bio.Seq.FASTQ
using Base.Test

import FASTQDemultiplexer: FASTQdemultiplex, Interpreter

# TODO replace the hand-written entries with a programatically
# generated strings.
# TODO test the readstring instead of comapring the data?


# we test the FASTQdemultiplexer on a mock-up protocol

cellbarcodes=DNASequence["AGTCCATGCT"]
umibarcodes=DNASequence["AAAAAT"]

interpreter = Interpreter(
    1,
    8:60,
    (4:7,1:6),
    (1:0,7:12),
    cellbarcodes,
    umibarcodes)

cellid = 1
umiid = 1
cellbarcode=cellbarcodes[cellid]
umibarcode=umibarcodes[umiid]

R1()=IOBuffer("""
@test:r1:1
AAACCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
123456789012345678901234567890123456789012345678901234567890
@test:r1:2
AAA$(cellbarcode[1:4])TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
123456789012345678901234567890123456789012345678901234567890
""")

R2()=IOBuffer("""
@test:r2:1
GGGGGGAAAAAA
+
123456789012
@test:r2:2
$(cellbarcode[5:10])$umibarcode
+
123456789012
""")


@testset "Basic functionality" begin

    cell1()=IOBuffer("""
    @test:r1:2
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +
    89012345678901234567890123456789012345678901234567890
    """)

    unmatched()=IOBuffer("""
    @test:r1:1
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +
    89012345678901234567890123456789012345678901234567890
    """)

    buffers=Dict()
    fastqstream(cellid) = buffers[cellid] = IOBuffer()

    FASTQdemultiplex((R1(),R2()),interpreter,celldescryptor=fastqstream,closebuffers=false)
    # rewind the buffers so we can read them again
    map(seekstart,values(buffers))

    @test haskey(buffers,cellid)
    @test haskey(buffers,0)

    @test buffers[0].data == unmatched().data
    @test buffers[1].data == cell1().data
end


@testset "Disabling the quality output" begin

    cell1()=IOBuffer("""
    @test:r1:2
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +

    """)

    unmatched()=IOBuffer("""
    @test:r1:1
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +

    """)

    buffers=Dict()
    fastqstream(cellid) = buffers[cellid] = IOBuffer()

    FASTQdemultiplex((R1(),R2()),interpreter,celldescryptor=fastqstream,closebuffers=false,outputquality=false)
    # rewind the buffers so we can read them again
    map(seekstart,values(buffers))

    @test haskey(buffers,cellid)
    @test haskey(buffers,0)

    @test buffers[0].data == unmatched().data
    @test buffers[1].data == cell1().data
end


@testset "UMI output" begin

    cell1()=IOBuffer("""
    @test:r1:2
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +
    89012345678901234567890123456789012345678901234567890
    """)

    unmatched()=IOBuffer("""
    @test:r1:1
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    +
    89012345678901234567890123456789012345678901234567890
    """)

    cell1umi()=IOBuffer("""
    $umibarcode
    """)

    unmatchedumi()=IOBuffer("""
    AAAAAA
    """)

    buffers=Dict()
    fastqstream(cellid) = buffers[cellid] = IOBuffer()

    umibuffers=Dict()
    umistream(cellid) = umibuffers[cellid] = IOBuffer()

    FASTQdemultiplex((R1(),R2()),interpreter,
                     celldescryptor=fastqstream,
                     umidescryptor=umistream,
                     closebuffers=false,
                     writeumis=true)
    # rewind the buffers so we can read them again
    map(seekstart,values(buffers))
    map(seekstart,values(umibuffers))

    @test haskey(buffers,cellid)
    @test haskey(buffers,0)

    @test buffers[0].data == unmatched().data
    @test buffers[1].data == cell1().data

    @test haskey(umibuffers,cellid)
    @test haskey(umibuffers,0)

    @test umibuffers[0].data == unmatchedumi().data
    @test umibuffers[1].data == cell1umi().data
end
