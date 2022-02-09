using RandomMatrixDistributions
using Random: seed!
using Test

@testset "RandomMatrixDistributions" begin
    seed!(1)

    @testset "MatrixSamplers" begin
        include("MatrixSamplers.jl")
    end
    
    @testset "Densities" begin
        include("Densities.jl")
    end
end
