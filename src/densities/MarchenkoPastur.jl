export MarchenkoPastur

"""
    MarchenkoPastur(γ::Real)

Marchenko-Pastur distribution, where 0 < `γ` ≤ 1.

The limiting spectral distribution of a p×p covariance matrix of n standard normal observations, where p/n → γ.
"""
struct MarchenkoPastur <: ContinuousUnivariateDistribution
    gamma::Real
    MarchenkoPastur(gamma) = 0 < gamma <= 1 ? new(gamma) : error("Gamma must be in (0, 1]")
end

function minimum(d::MarchenkoPastur)
    (1 - sqrt(d.gamma))^2
end

function maximum(d::MarchenkoPastur)
    (1 + sqrt(d.gamma))^2
end

function pdf(d::MarchenkoPastur, x::Real)
    lambdamin = minimum(d)
    lambdamax = maximum(d)

    if lambdamin < x < lambdamax
        return sqrt((lambdamax - x) * (x - lambdamin))/(d.gamma * x * 2 * pi)
    else
        return 0
    end
end
