using RandomMatrixDistributions
using Distributions
using Test

@testset "WhiteSamplers" begin

    n = 1000
    p = 400

    # test white Wishart
    for beta in [1, 4, 8]
        #println("Testing white Wishart with β = ", beta)
        
        λs = randeigvals(SpikedWishart(beta, n, p))

        @test length(λs) == p
        @test eltype(λs) <: Real

        @test minimum(λs) > 0

        #println("Testing Jacobi with β = ", beta)
        λs = randeigvals(Jacobi(beta, n, n, p))

        @test length(λs) == p
        @test eltype(λs) <: Real

        @test minimum(λs) > 0
    end

end

@testset "TridiagSamplers" begin

    n = 1000
    p = 400 

    for beta in [1, 2]
        #println("Testing Spiked Wishart with β = 1")
        s = 8
        
        for W in [SpikedWishart(beta, n, p, [s], scaled = false), SpikedWigner(beta, p, [s], scaled = false)]
            λs = randeigvals(W)
            
            @test length(λs) == p
            @test eltype(λs) <: Real

            # Check if the maximum spike is reasonable
            # This tests can fail in principle, but only with negligible probability
            spikedist = supercrit_dist(W)
            μ = mean(spikedist)[end]
            σ = sqrt(var(spikedist)[end])

            @test abs(maximum(λs) - μ)/σ < 4

            @test all(randeigstat(W, eigmax, 10) .> 0)
        end
    end

end

@testset "BandedSamplers" begin

    n = 1000
    p = 400 

    for beta in [1, 2]
        #println("Testing Spiked Wishart with β = 1")
        s = 8
        for W in [SpikedWishart(beta, n, p, [0.1, 2, 4, s], scaled = true), SpikedWigner(beta, p, [0.1, 2, 4, s], scaled = true)]
            λs = randeigvals(W)
            
            @test length(λs) == p
            @test eltype(λs) <: Real

            # Check if the maximum spike is reasonable
            # This tests can fail in principle, but only with negligible probability
            spikedist = supercrit_dist(W)
            μ = mean(spikedist)[end]
            σ = sqrt(var(spikedist)[end])

            @test abs(maximum(λs) - μ)/σ < 4

            @test all(randeigstat(W, eigmax, 10) .> 0)
        end
    end

end
