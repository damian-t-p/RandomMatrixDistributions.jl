export TracyWidom

struct TracyWidom <: ContinuousUnivariateDistribution
    beta::Integer
end

_TWfs = Dict()

function cdf(d::TracyWidom, x)
    cdf(RandomMatrices.TracyWidom, x, beta=d.beta)
end

function _TWf(x, beta)
    if !(beta in keys(_TWfs))
        dom = Chebyshev(-5..5)
        D = Derivative(dom)
        _TWfs[beta] = D * Fun(x -> cdf(RandomMatrices.TracyWidom, x, beta=beta), dom)
    end
    _TWfs[beta](x)
end

function pdf(d::TracyWidom, x)
    _TWf(x, d.beta)
end

