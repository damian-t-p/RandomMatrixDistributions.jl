using RandomMatrixDistributions
using Test

n = 1000
p = 400
γ = p/n

# test white Wishart
for beta in [1, 4, 8]
    λs = randeigvals(SpikedWishart(beta, n, p))

    @test length(λs) == p
    @test eltype(λs) <: Real

    @test minimum(λs) > 0

    λs = randeigvals(Jacobi(beta, n, n, p))

    @test length(λs) == p
    @test eltype(λs) <: Real

    @test minimum(λs) > 0
end

s = 8
λs = randeigvals(SpikedWishart(1, n, p, [2, 4, s]))/n

@test length(λs) == p
@test eltype(λs) <: Real

# Check if the maximum spike is reasonable
# This tests can fail in principle, but only with negligible probability
μ = (1 + s) * (1 + γ/s)
σ = (1 + s) * sqrt(γ * 2 * (1 - γ/s^2))

@test abs(maximum(λs) - μ)/σ < 4
