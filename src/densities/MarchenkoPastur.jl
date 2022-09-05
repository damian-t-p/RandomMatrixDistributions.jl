export MarchenkoPastur

"""
    MarchenkoPastur(γ::Real)

Marchenko-Pastur distribution, where 0 < `γ`.

The limiting spectral distribution of a p×p covariance matrix of n standard normal observations, where p/n → γ.
"""
struct MarchenkoPastur <: ContinuousUnivariateDistribution
    gamma::Real
    MarchenkoPastur(gamma) = 0 < gamma ? new(gamma) : error("Gamma must be > 0")
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
    elseif x == 0.0 && d.gamma > 1.0
        1. - 1.0/d.gamma
    else
        return 0.0
    end
end
