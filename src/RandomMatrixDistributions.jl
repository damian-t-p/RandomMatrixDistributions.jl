module RandomMatrixDistributions

import RandomMatrices

using Random, Distributions
using LinearAlgebra
using BandedMatrices
using KrylovKit
using ApproxFun
using Interpolations

import LinearAlgebra: eigmax

import Distributions: minimum, maximum, quantile

import Distributions: ContinuousUnivariateDistribution,
    ContinuousMatrixDistribution,
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment


export randeigvals, randeigstat,
    minimum, maximum, quantile
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment,
    eigmax


# Generic eigenvalue sampling functions

function randeigvals(d)
    eigvals(randreduced(d))
end

function randeigstat(d, eigstat, nsims::Int)
    statvals = Array{Float64, 1}(undef, nsims)
    for i in 1:nsims
        statvals[i] = eigstat(randreduced(d))
    end
    statvals
end

include("BandedHelpers.jl")

include("SpikedWishart.jl")

include("Jacobi.jl")

include("densities/MarchenkoPastur.jl")
include("densities/TracyWidom.jl")
include("densities/Wachter.jl")

end
