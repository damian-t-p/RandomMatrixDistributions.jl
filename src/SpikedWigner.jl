export SpikedWigner

"""
    SpikedWigner(β::Int, n::Int[, spikes::Vector{Float64}, scaled::Bool=false])

Distribution on an n×n spiked Gaussian Wigner matrix.

Wigner matrices are Hermitian with independent real, complex or quaternion standard Gaussian entries depending on whether β = 1, 2 or 4.

A diagonal matrix with entries given by `spikes` multiplied by `√n` is added to produce the spiked Wigner matrix.

If `scaled == true`, then the resulting matrix is divided by `√n` so that its bulk distribution converges to the semicircle law supported on [-2, 2].
"""
struct SpikedWigner <: Distribution{Matrixvariate, ComplexContinuous}
    beta::Integer
    n::Integer
    spikes::Vector{Float64}
    scaled::Bool
end

SpikedWigner(beta, n, spikes; scaled=false) = SpikedWigner(beta, n, spikes, scaled)
SpikedWigner(beta, n; scaled=false) = SpikedWigner(beta, n, [], scaled)

# PROPERTIES

function scaling(d::SpikedWigner)
    d.scaled ? 1/sqrt(d.n) : 1
end

Base.size(d::SpikedWigner) = (d.n, d.n)

function eltype(d::SpikedWigner)
    if d.beta == 1
        Float64
    elseif d.beta == 2
        ComplexF64
    else
        error("Wigner matrices only instantiated for β = 1, 2")
    end
end

# SAMPLERS

function _rand!(rng::AbstractRNG, d::SpikedWigner, x::DenseMatrix{T}) where {T <: Number}
    if d.beta == 1
        A = randn(rng, d.n, d.n)
    elseif d.beta == 2
        A = (randn(rng, d.n, d.n) + im * randn(rng, d.n, d.n))/sqrt(2)
    else
        error("Wigner matrices only instantiated for β = 1, 2")
    end

    x[:] = (A + A')/sqrt(2)

    for (i, spike) in enumerate(d.spikes)
        x[i,i] += sqrt(d.n) * spike
    end
    
    x *= scaling(d)

    return x
    
end

function randreduced(rng::AbstractRNG, d::SpikedWigner)
    if length(d.spikes) <= 1
        return randtridiagonal(rng, d)
    else
        return randbanded(rng, d)
    end
end

function randtridiagonal(rng::AbstractRNG, d::SpikedWigner)
    
    # diagonal and superdiagonal of B
    dv = rand(rng, Normal(0, sqrt(2)), d.n)/sqrt(d.beta)
    ev = [rand(rng, Chi(d.beta*(d.n-k))) for k in 1:(d.n-1)]/sqrt(d.beta)

    if length(d.spikes) == 1
        dv[1] += d.spikes[1] * sqrt(d.n)
    end

    SymTridiagonal(dv, ev) * scaling(d)
end

function randbanded(rng::AbstractRNG, d::SpikedWigner)
    @assert d.beta in [1, 2]
    r = length(d.spikes)

    if d.beta == 1
        U = BandedMatrix{Float64}(undef, (d.n, d.n), (r, r))
    elseif d.beta == 2
        U = BandedMatrix{Complex{Float64}}(undef, (d.n, d.n), (r, r))
    end

    dv = randn(rng, d.n) * sqrt(2/d.beta)
    @. dv[1:r] += d.spikes * sqrt(d.n)

    U[band(0)] .= dv 
    
    for k = 1:(r-1)
        if d.beta == 1
            ev = randn(rng, d.n - k)
        elseif d.beta == 2
            ev = (randn(rng, d.n - k) + im * randn(rng, d.n-k))/sqrt(2)
        end

        U[band(k)] .= ev
    end

    U[band(r)] .= [rand(rng, Chi(d.beta*(d.n - k)))/sqrt(d.beta) for k in r:(d.n-1)]

    if d.beta == 1
        Symmetric(U) * scaling(d)
    elseif d.beta == 2
        # The conjugate transpose is done like this rather than with ' because
        # U'U is not automatically a banded matrix
        Hermitian(U) * scaling(d)
    end
        
end

function supercrit_spikes(d::SpikedWigner)
    d.spikes[d.spikes .> 1]
end

# For beta = 1, see Feral & Peche 2006 Theorem 1.1
# For beta = 2, see Peche Theorem 1.2
function supercrit_dist(d::SpikedWigner)

    cspikes = supercrit_spikes(d)
    
    if !allunique(cspikes)
        throw("Supercritical spikes with multiplicity > 1 not supported")
    end

    mu = @. cspikes + 1/cspikes
    sigma = @. sqrt(2/d.beta) * sqrt(cspikes^2 - 1) / (cspikes * sqrt(d.n))

    if d.scaled == false
        mu *= sqrt(d.n)
        sigma *= sqrt(d.n)
    end
    
    MvNormal(mu, Diagonal(sigma .^ 2))
end

