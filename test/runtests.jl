using RandomMatrixDistributions
using Random: seed!
using Test

@testset "RandomMatrixDistributions" begin
    seed!(1)
    include("MatrixSamplers.jl")
end
