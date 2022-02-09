using RandomMatrixDistributions
using Distributions
using Random
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
        for W in [
            SpikedWishart(beta, n, p, [0.1, 2, 4, s], scaled = true),
            SpikedWigner(beta, p, [0.1, 2, 4, s], scaled = true)]

            λs = rand(EigvalDist(W))
            
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

@testset "DenseSamplers" begin
    
    n = 100
    p = 40 

    for beta in [1, 2, 4]
        #println("Testing Spiked Wishart with β = 1")
        s = 8
        for W in [
            SpikedWishart(beta, n, p, [0.1, 2, 4, s], scaled = true),
            SpikedWigner(beta, p, [0.1, 2, 4, s], scaled = true),
            Jacobi(beta, n, n, p)]

            if beta == 4
                @test_throws ErrorException rand(W)
            else
                M = rand(W)

                if beta == 1
                    @test eltype(M) <: Real
                elseif beta == 2
                    @test eltype(M) <: Complex
                end
                
            end
        end
    end
    
end

@testset "RandomnessReproducibility" begin

    n = 10
    p = 5

    beta = 1
    
    # white
    for W in [
        SpikedWishart(beta, n, p, scaled = true),
        SpikedWigner(beta, p, scaled = true),
        Jacobi(beta, n, n, p)]

        λ1 = rand(MersenneTwister(0), EigvalDist(W))
        λ2 = rand(MersenneTwister(0), EigvalDist(W))
        
        @test λ1 ≈ λ2

        M1 = rand(MersenneTwister(0), W)
        M2 = rand(MersenneTwister(0), W)

        @test M1 ≈ M2
    end
    
    # tridiagonal
    s = 8
    for W in [
        SpikedWishart(beta, n, p, [s], scaled = true),
        SpikedWigner(beta, p, [s], scaled = true)]

        λ1 = rand(MersenneTwister(0), EigvalDist(W))
        λ2 = rand(MersenneTwister(0), EigvalDist(W))
        
        @test λ1 ≈ λ2

        M1 = rand(MersenneTwister(0), W)
        M2 = rand(MersenneTwister(0), W)

        @test M1 ≈ M2
    end
    
    # banded
    s = 8
    for W in [
        SpikedWishart(beta, n, p, [0.1, 2, 4, s], scaled = true),
        SpikedWigner(beta, p, [0.1, 2, 4, s], scaled = true)]

        λ1 = rand(MersenneTwister(0), EigvalDist(W))
        λ2 = rand(MersenneTwister(0), EigvalDist(W))
        
        @test λ1 ≈ λ2

        M1 = rand(MersenneTwister(0), W)
        M2 = rand(MersenneTwister(0), W)

        @test M1 ≈ M2
    end
    
end

