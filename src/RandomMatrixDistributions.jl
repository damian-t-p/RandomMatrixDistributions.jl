module RandomMatrixDistributions

using Random, Distributions
using LinearAlgebra
using BandedMatrices
using KrylovKit

import LinearAlgebra: eigmax

import Distributions: minimum, maximum, quantile

import Distributions: ContinuousUnivariateDistribution,
    ContinuousMatrixDistribution,
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment


export randeigvals, randeigstat, supercrit_dist
    minimum, maximum, quantile
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment,
    eigmax


# Generic eigenvalue sampling functions

"""
    randeigvals(d::ContinuousMatrixDistribution)

Sample a vector of eigenvalues of a matrix drawn from the matrix ensemble `d`.
"""
function randeigvals(d::ContinuousMatrixDistribution)
    eigvals(randreduced(d))
end

"""
    randeigstat(d::ContinuousMatrixDistribution, eigstat::Function, n::Int)

Sample `n` realisations of the eigenvalue statistic `eigstat` evaluated at a matrices drawn from the ensemble `d`.

`eigstat` is a function of a square matrix argument whose value depends only on the eigenvalues of that matrix.

# Usage

```julia
rangeigstat(SpikedWigner(2, 50), eigmax, 100)
```
"""
function randeigstat(d::ContinuousMatrixDistribution, eigstat::Function, n::Int)
    statvals = Array{Float64, 1}(undef, n)
    for i in 1:n
        statvals[i] = eigstat(randreduced(d))
    end
    statvals
end


"""
    supercrit_dist(d::ContinuousMatrixDistribution)

Compute the approximate joint distribution of the supercritical eigenvalues of the ensemble `d`.
"""
function supercrit_dist(d::ContinuousMatrixDistribution)
end

include("BandedHelpers.jl")

include("SpikedWigner.jl")
include("SpikedWishart.jl")

include("Jacobi.jl")

include("densities/MarchenkoPastur.jl")
include("densities/TracyWidom.jl")
include("densities/Wachter.jl")

end
