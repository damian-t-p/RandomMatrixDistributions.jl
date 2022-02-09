export SpikedWishart

"""
    SpikedWishart(β::Int, n::Int, p::Int[, spikes::Vector{Float64}, scaled::Bool=false)

Distribution of a p×p spiked Wishart matrix.

If X is a p×n matrix with independent real, complex or quaternion standard Gaussian entries, depending on whether β = 1, 2 or 4, then XX† has a `Wishart(β, n, p)` distribution.

If Λ is a diagonal matrix whose entries are `√(1 .+ spikes)`, then ΛXX†Λ has a `SpikedWishart(β, n, p, spikes)` distribution.

If `scaled == true`, then the resulting matrix is divided by `p` so that its bulk distribution converges to the Marchenko-Pastur law.
"""
struct SpikedWishart <: Distribution{Matrixvariate, ComplexContinuous}
    beta::Integer
    n::Integer
    p::Integer
    spikes::Vector{Float64}
    scaled::Bool
end

SpikedWishart(beta, n, p, spikes; scaled=false) = SpikedWishart(beta, n, p, spikes, scaled)
SpikedWishart(beta, n, p; scaled=false) = SpikedWishart(beta, n, p, [], scaled)

# PROPERTIES

function scaling(d::SpikedWishart)
    d.scaled ? 1/d.n : 1
end

Base.size(d::SpikedWishart) = (d.p, d.p)

function eltype(d::SpikedWishart)
    if d.beta == 1
        Float64
    elseif d.beta == 2
        ComplexF64
    else
        error("Wishart matrices only instantiated for β = 1, 2")
    end
end

# SAMPLERS

function _rand!(rng::AbstractRNG, d::SpikedWishart, x::DenseMatrix{T}) where {T <: Number}
    if d.beta == 1
        A = randn(rng, d.p, d.n)
    elseif d.beta == 2
        A = (randn(rng, d.p, d.n) + im * randn(rng, d.p, d.n))/sqrt(2)
    else
        error("Wishart matrices only instantiated for β = 1, 2")
    end

    for (i, spike) in enumerate(d.spikes)
        A[i,:] *= sqrt(1 + spike)
    end
    
    x[:] = (A * A') * scaling(d)
    
    return x
    
end

function randreduced(rng::AbstractRNG, d::SpikedWishart)
    if length(d.spikes) <= 1
        return randtridiagonal(rng, d)
    else
        return randbanded(rng, d)
    end
end

function randtridiagonal(rng::AbstractRNG, d::SpikedWishart)
    a = d.beta * d.n/2
    
    # diagonal and superdiagonal of B
    Bdv = [rand(rng, Chi(2*a - d.beta*k))/sqrt(d.beta) for k in 0:(d.p-1)]
    Bev = [rand(rng, Chi(d.beta*(d.p-k)))/sqrt(d.beta) for k in 1:(d.p-1)]

    if length(d.spikes) == 1
        Bdv[1] *= sqrt(1 + d.spikes[1])
    end
    
    # diagonal and off-diagonals of BB'
    dv = Bdv.^2 + [0; Bev.^2]
    ev = Bdv[1:end-1] .* Bev
    
    SymTridiagonal(dv, ev) * scaling(d)
end

# Multi-spike version only implemented for beta=1
function randbanded(rng::AbstractRNG, d::SpikedWishart)
    @assert d.beta in [1, 2]
    r = length(d.spikes)

    if d.beta == 1
        U = BandedMatrix{Float64}(undef, (d.p, d.p), (0, r))
    elseif d.beta == 2
        U = BandedMatrix{Complex{Float64}}(undef, (d.p, d.p), (0, r))
    end

    dv = [rand(rng, Chi(d.beta*(d.n - k + 1)))/sqrt(d.beta) for k in 1:d.p]
    @. dv[1:r] *= sqrt(1 + d.spikes)

    U[band(0)] .= dv
    
    for k = 1:(r-1)
        if d.beta == 1
            ev = randn(rng, d.p - k)
        elseif d.beta == 2
            ev = (randn(rng, d.p - k) + im * randn(rng, d.p-k))/sqrt(2)
        end
        
        @. ev[1:(r-k)] *= sqrt(1 + d.spikes[(k+1):end])

        U[band(k)] .= ev
    end

    U[band(r)] .= [rand(rng, Chi(d.beta*(d.p - k)))/sqrt(d.beta) for k in r:(d.p-1)]

    if d.beta == 1
        Symmetric(U' * U) * scaling(d)
    elseif d.beta == 2
        # The conjugate transpose is done like this rather than with ' because
        # U'U is not automatically a banded matrix
        Hermitian(transpose(conj(U)) * U) * scaling(d)
    end
        
end

function supercrit_spikes(d::SpikedWishart)
    gamma = d.p/d.n
    d.spikes[d.spikes .> sqrt(gamma)]
end


# For beta = 1, see Paul 2007 Theorem 3
# For beta = 2, see BBP 2005 Theorem 1.1(b)
function supercrit_dist(d::SpikedWishart)

    gamma = d.p/d.n
    cspikes = supercrit_spikes(d)

    if !allunique(cspikes)
        throw("Supercritical spikes with multiplicity > 1 not supported")
    end
    
    l = 1 .+ cspikes
    
    mu = @. l * (1 + gamma/cspikes)
    sigma = @. sqrt(2/d.beta) * l * sqrt(gamma  * (1 - gamma/cspikes^2)) / sqrt(d.n)

    if d.scaled == false
        mu *= d.n
        sigma *= d.n
    end
    
    MvNormal(mu, Diagonal(sigma .^ 2))
end
