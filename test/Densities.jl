using RandomMatrixDistributions
using Distributions
using Test

@testset "LimitingSpectralDistributions" begin

    dist = SpikedWigner(1, 10, [10])
    @test bulk_dist(dist) isa Distributions.Semicircle
    @test bulk_dist(dist) == bulk_dist(EigvalDist(dist))
    @test supercrit_dist(dist) == supercrit_dist(EigvalDist(dist))
    
    @test bulk_dist(SpikedWishart(2, 10, 5)) isa MarchenkoPastur
    @test bulk_dist(Jacobi(2, 10, 10, 5)) isa Wachter
    
end

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

@testset "TracyWidom" begin

    @test_throws ErrorException TracyWidom(10)

    # computed from the Mathematica function TracyWidomDistribution
    cdf_vals = Dict(1 => 0.584, 2 => 0.807, 4 => 0.961)
    pdf_vals = Dict(1 => 0.304, 2 => 0.286, 4 => 0.107)
    qnt_vals = Dict(1 => -1.269, 2 => -1.805, 4 => -2.327)
    
    for beta in [1, 2, 4]

        dist = TracyWidom(beta)

        @test isapprox(cdf(dist, -1), cdf_vals[beta], atol = 1e-3)
        @test isapprox(pdf(dist, -1), pdf_vals[beta], atol = 1e-3)
        @test isapprox(quantile(dist, 0.5), qnt_vals[beta], atol = 1e-3)

    end
    
end
