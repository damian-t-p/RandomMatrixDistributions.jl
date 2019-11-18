using RandomMatrixDistributions
using Test

@testset "MarchenkoPastur" begin

    @test_throws ErrorException MarchenkoPastur(5)

    dist = MarchenkoPastur(1)

    @test minimum(dist) == 0
    @test maximum(dist) == 4
    
    @test pdf(dist, -1) == pdf(dist, 5) == 0
    @test pdf(dist, 2) ≈ 1/(2*pi)
end

@testset "Wachter" begin

    @test_throws ErrorException Wachter(10, 10)

    dist = Wachter(1, 0.5)

    @test minimum(dist) == 0
    @test maximum(dist) == 16

    @test pdf(dist, -1) == pdf(dist, 17) == 0
    @test pdf(Wachter(0.5, 0), 1) ≈ pdf(MarchenkoPastur(0.5), 1)
end

