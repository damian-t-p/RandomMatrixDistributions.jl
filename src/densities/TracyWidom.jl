export TracyWidom

struct TracyWidom <: ContinuousUnivariateDistribution
    beta::Integer
end

function cdf(d::TracyWidom, x::Real)
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

function pdf(d::TracyWidom, x::Real)
    _TWf(x, d.beta)
end

_TWinvs = Dict()

function _TWinv(t, beta)
    (0 < t < 1) || throw(DomainError(t))

    if !(beta in keys(_TWinvs))
        # bounds outside which we use limiting representation
        u = 4
        l = -5
        
        # TODO: make a more principled grid using tail behaviour
        n = 1000
        y = range(l, u, length=n)
        x = cdf(RandomMatrices.TracyWidom, y, beta=1)

        # Transform limiting ends to maintain differentiability
        righta = pdf(TracyWidom(beta), u)/(beta * sqrt(u) * exp(-beta * u^(3/2)*2/3))
        rightb = x[end] - righta * (1 - exp(-beta * u^(3/2) * 2/3))

        lefta = pdf(TracyWidom(beta), l)/(beta * l^2/8 * exp(beta * l^3/24))
        leftb = x[1] - lefta * exp(-beta * abs(l)^3/24)

        itp = interpolate((x,), y, Gridded(Linear()))

        _TWinvs[beta] = [righta, rightb, x[end], lefta, leftb, x[1], itp]
    end

    # TODO: make the indexing less terrible
    if t > _TWinvs[beta][3]
        s = (t - _TWinvs[beta][2])/_TWinvs[beta][1]
        return (log(1/(1 - s)) * 3/(2 * beta))^(2/3)
    elseif t < _TWinvs[beta][6]
        s = (t - _TWinvs[beta][5])/_TWinvs[beta][4]
        return -(log(1/s) * 24/beta)^(1/3)
    else
        return _TWinvs[beta][end](t)
    end
        
end

function quantile(d::TracyWidom, q::Real)
    _TWinv(q, d.beta)
end

