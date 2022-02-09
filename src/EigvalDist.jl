export EigvalDist

"""
    EigvalDist(matrixdist::ContinuousMatrixDistribution)

The joint distribution of eigenvalues sampled from a random matrix with the distribution `matrixdist`.
"""
struct EigvalDist <: ContinuousMultivariateDistribution
    matrixdist::ContinuousMatrixDistribution
end

Base.length(d::EigvalDist) = Base.size(d.matrixdist)[1]

function _rand!(rng::AbstractRNG, d::EigvalDist, x::AbstractVector{T}) where {T <: Real}
    x[:] = randeigvals(rng, d.matrixdist)
    return x
end
