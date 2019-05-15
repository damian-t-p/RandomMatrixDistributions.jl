export MarchenkoPastur

struct MarchenkoPastur <: ContinuousUnivariateDistribution
    gamma::Real
    MarchenkoPastur(gamma) = 0 < gamma < 1 ? new(gamma) : error("Gamma must be in (0, 1)")
end

function minimum(d::MarchenkoPastur)
    (1 - sqrt(d.gamma))^2
end

function maximum(d::MarchenkoPastur)
    (1 + sqrt(d.gamma))^2
end

function pdf(d::MarchenkoPastur, x)
    lambdamin = minimum(d)
    lambdamax = maximum(d)

    if lambdamin < x < lambdamax
        return sqrt((lambdamax - x) * (x - lambdamin))/(d.gamma * x * 2 * pi)
    else
        return 0
    end
end
