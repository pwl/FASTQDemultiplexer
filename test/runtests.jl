using FASTQDemultiplexer
using Base.Test

import FASTQDemultiplexer: idtoname, gen_id, Barcodes, getgroup

# utils
@testset "utils" begin
    ain = UInt8['A','C','T','A','C','A','T','G','A','C']
    @test idtoname(gen_id(ain),length(ain))==String(ain)
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

    barcodes = Barcodes([group1,group2],["group1","group2"])

    @test getgroup(barcodes,gen_id(ain1)) == (1, "group1")
    @test getgroup(barcodes,gen_id(ain4)) == (2, "group2")
    @test getgroup(barcodes,gen_id(ainnone)) == (-1, "unassigned")
end
