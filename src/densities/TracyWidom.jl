export TracyWidom

struct TracyWidom <: ContinuousUnivariateDistribution
    beta::Integer
end

function cdf(d::TracyWidom, x)
    cdf(RandomMatrices.TracyWidom, x, beta=d.beta)
end

_TWfs = Dict()

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

_TWinvs = Dict()

function _TWinv(t, beta)
    if !(beta in keys(_TWinvs))
        # TODO: make a more principled grid using tail behaviour
        n = 1000
        y = range(-5, 4, length=n)
        x = cdf(RandomMatrices.TracyWidom, y, beta=1)

        itp = interpolate((x,), y, Gridded(Linear()))
        _TWinvs[beta] = Interpolations.extrapolate(itp, Interpolations.Line())
    end
    _TWinvs[beta](t)
end

function quantile(d::TracyWidom, q::Real)
    _TWinv(q, d.beta)
end

