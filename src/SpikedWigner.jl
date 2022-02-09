export SpikedWigner

"""
    SpikedWigner(β::Int, n::Int[, spikes::Vector{Float64}, scaled::Bool=false])

Distribution on an n×n spiked Gaussian Wigner matrix.

Wigner matrices are Hermitian with independent real, complex or quaternion standard Gaussian entries depending on whether β = 1, 2 or 4.

A diagonal matrix with entries given by `spikes` multiplied by `√n` is added to produce the spiked Wigner matrix.

If `scaled == true`, then the resulting matrix is divided by `√n` so that its bulk distribution converges to the semicircle law supported on [-2, 2].
"""
struct SpikedWigner <: ContinuousMatrixDistribution
    beta::Integer
    n::Integer
    spikes::Vector{Float64}
    scaled::Bool
end

SpikedWigner(beta, n, spikes; scaled=false) = SpikedWigner(beta, n, spikes, scaled)
SpikedWigner(beta, n; scaled=false) = SpikedWigner(beta, n, [], scaled)

# SAMPLERS

function scaling(d::SpikedWigner)
    d.scaled ? 1/sqrt(d.n) : 1
end

function randreduced(d::SpikedWigner)
    if length(d.spikes) <= 1
        return randtridiagonal(d)
    else
        return randbanded(d)
    end
end

function randtridiagonal(d::SpikedWigner)
    
    # diagonal and superdiagonal of B
    dv = rand(Normal(0, sqrt(2)), d.n)/sqrt(d.beta)
    ev = [rand(Chi(d.beta*(d.n-k))) for k in 1:(d.n-1)]/sqrt(d.beta)

    if length(d.spikes) == 1
        dv[1] += d.spikes[1] * sqrt(d.n)
    end

    SymTridiagonal(dv, ev) * scaling(d)
end

function randbanded(d::SpikedWigner)
    @assert d.beta in [1, 2]
    r = length(d.spikes)

    if d.beta == 1
        U = BandedMatrix{Float64}(undef, (d.n, d.n), (r, r))
    elseif d.beta == 2
        U = BandedMatrix{Complex{Float64}}(undef, (d.n, d.n), (r, r))
    end

    dv = randn(d.n) * sqrt(2/d.beta)
    @. dv[1:r] += d.spikes * sqrt(d.n)

    U[band(0)] .= dv 
    
    for k = 1:(r-1)
        if d.beta == 1
            ev = randn(d.n - k)
        elseif d.beta == 2
            ev = (randn(d.n - k) + im * randn(d.n-k))/sqrt(2)
        end

        U[band(k)] .= ev
    end

    U[band(r)] .= [rand(Chi(d.beta*(d.n - k)))/sqrt(d.beta) for k in r:(d.n-1)]

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

