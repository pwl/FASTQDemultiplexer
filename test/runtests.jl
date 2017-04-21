using FASTQDemultiplexer
using Bio.Seq
using Bio.Seq.FASTQ
using Base.Test

import FASTQDemultiplexer: iterate, Interpreter


# we test the demultiplexer on the MARSSeq protocol

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

R1()=IOBuffer(
"""
@test:r1:1
AAACCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
123456789012345678901234567890123456789012345678901234567890
@test:r1:2
AAA$(cellbarcode[1:4])TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
123456789012345678901234567890123456789012345678901234567890
""")

R2()=IOBuffer(
"""
@test:r2:1
GGGGGGAAAAAA
+
123456789012
@test:r2:2
$(cellbarcode[5:10])$umibarcode
+
123456789012
""")

cell1()=IOBuffer(
"""
@test:r1:2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
89012345678901234567890123456789012345678901234567890
""")

unmatched()=IOBuffer(
"""
@test:r1:1
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
89012345678901234567890123456789012345678901234567890
""")

cell1_noqual()=IOBuffer(
"""
@test:r1:2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+

""")

unmatched_noqual()=IOBuffer(
"""
@test:r1:1
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+

""")


buffers=Dict()
fastqstream(cellid) = buffers[cellid] = IOBuffer()

iterate((R1(),R2()),interpreter,gendescryptor=fastqstream,closebuffers=false)
# rewind the buffers so we can read them again
map(seekstart,values(buffers))

@test haskey(buffers,cellid)
@test haskey(buffers,0)

@test buffers[0].data == unmatched().data
@test buffers[1].data == cell1().data

buffers=Dict()
fastqstream(cellid) = buffers[cellid] = IOBuffer()

iterate((R1(),R2()),interpreter,gendescryptor=fastqstream,closebuffers=false,outputquality=false)
# rewind the buffers so we can read them again
map(seekstart,values(buffers))

@test haskey(buffers,cellid)
@test haskey(buffers,0)

@test buffers[0].data == unmatched_noqual().data
@test buffers[1].data == cell1_noqual().data
