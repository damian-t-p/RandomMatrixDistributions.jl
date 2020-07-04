using FastGaussQuadrature
using SpecialFunctions
using ApproxFun
using Interpolations

export TracyWidom

struct TracyWidom <: ContinuousUnivariateDistribution
    beta::Integer
end


function cdf(d::TracyWidom, x::Real)
    if d.beta == 1
        _fredholm_det(_V_airy_kernel, x, 100)
    elseif d.beta == 2
        _fredholm_det(_airy_kernel, x, 100)
    elseif d.beta == 4
        (_fredholm_det(_V_airy_kernel, x * sqrt(2), 100) + _fredholm_det(_V_airy_kernel, x * sqrt(2), 100, -1))/2
    else
        throw(DomainError(d.beta))
    end
end

_TWfs = Dict()

function _TWf(x, beta)
    if !(beta in keys(_TWfs))
        dom = Chebyshev(-5..5)
        D = Derivative(dom)
        _TWfs[beta] = D * Fun(x -> cdf(TracyWidom(beta), x), dom)
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
        x = cdf.(TracyWidom(beta), y)

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

function _phi(s, xi)
    s + 10 * tan(pi * (xi + 1)/4)
end

function _Dphi(s, xi)
    pi * 5/2 * sec(pi * (xi + 1)/4)^2
end

"""
Transform a kernel on L_2(s, infinity) to one on L_2(-1, 1)
"""
function _transform_kernel(kernel, s)
    (x, y) -> sqrt(_Dphi(s, x) * _Dphi(s, y)) * kernel(_phi(s, x), _phi(s, y))
end

"""
Approximate fredholm determinant of a kernel on (s, infinity)
"""
function _fredholm_det(kernel, s, N, z=1)
    
    (nodes, weights) = gausslegendre(N)
    trans_kernel = _transform_kernel(kernel, s)

    M = I - z * sqrt.(weights * weights') .* trans_kernel.(nodes', nodes)
    det(M)
end

function _airy_kernel(x, y)
    if x == y
        airyaiprime(x)^2 - x * airyai(x)^2
    else
        (airyai(x) * airyaiprime(y) - airyaiprime(x) * airyai(y))/(x - y)
    end
end

function _V_airy_kernel(x, y)
    airyai((x + y)/2)/2
end


