module RandomMatrixDistributions

using Random, Distributions
using LinearAlgebra
using BandedMatrices
using KrylovKit

import LinearAlgebra: eigmax

import Base: length, size, eltype

import Distributions: ContinuousUnivariateDistribution,
    ContinuousMatrixDistribution,
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment,
    minimum, maximum, quantile,
    _rand!

import Random: rand

export randeigvals, randeigstat, supercrit_dist, bulk_dist
    minimum, maximum, quantile,
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment,
    eigmax


struct ComplexContinuous <: ValueSupport end

include("EigvalDist.jl")

# Generic eigenvalue sampling functions

"""
    randeigvals([rng::AbstractRNG, ]d::MatrixDistribution)

Sample a vector of eigenvalues of a matrix drawn from the matrix ensemble `d`.
"""
function randeigvals(rng::AbstractRNG, d::MatrixDistribution)
    eigvals(randreduced(rng, d))
end

randeigvals(d::MatrixDistribution) = randeigvals(Random.default_rng(), d)

"""
    randeigstat([rng::AbstractRNG, ]d::MatrixDistribution, eigstat::Function, n::Int)

Sample `n` realisations of the eigenvalue statistic `eigstat` evaluated at a matrices drawn from the ensemble `d`.

`eigstat` is a function of a square matrix argument whose value depends only on the eigenvalues of that matrix.

# Usage

```julia
rangeigstat(SpikedWigner(2, 50), eigmax, 100)
```
"""
function randeigstat(rng::AbstractRNG, d::MatrixDistribution, eigstat::Function, n::Int)
    statvals = Array{Float64, 1}(undef, n)
    for i in 1:n
        statvals[i] = eigstat(randreduced(rng, d))
    end
    statvals
end

randeigstat(d::MatrixDistribution, eigstat::Function, n::Int) = randeigstat(Random.default_rng(), d, eigstat, n)

"""
    bulk_dist(d::Union{MatrixDistribution, EigvalDist})

Compute the limiting spectral distribution of `d`.
"""
bulk_dist(d::MatrixDistribution)
bulk_dist(d::EigvalDist) = bulk_dist(d.matrixdist)

"""
    supercrit_dist(d::Union{MatrixDistribution, EigvalDist})

Compute the approximate joint distribution of the supercritical eigenvalues of the ensemble `d`.
"""
supercrit_dist(d::MatrixDistribution)
supercrit_dist(d::EigvalDist) = supercrit_dist(d.matrixdist)

include("BandedHelpers.jl")

include("SpikedWigner.jl")
include("SpikedWishart.jl")
include("Jacobi.jl")

include("densities/MarchenkoPastur.jl")
include("densities/TracyWidom.jl")
include("densities/Wachter.jl")

end
