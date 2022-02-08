export SpikedWishart, supercrit_dist

struct SpikedWishart <: ContinuousMatrixDistribution
    beta::Integer
    n::Integer
    p::Integer
    spikes::Array{Float64, 1}
    scaled::Bool
end

SpikedWishart(beta, n, p, spikes; scaled=false) = SpikedWishart(beta, n, p, spikes, scaled)
SpikedWishart(beta, n, p; scaled=false) = SpikedWishart(beta, n, p, [], scaled)

# SAMPLERS

function scaling(d::SpikedWishart)
    d.scaled ? 1/d.n : 1
end

function randreduced(d::SpikedWishart)
    if length(d.spikes) <= 1
        return randtridiagonal(d)
    else
        return randbanded(d)
    end
end

function randtridiagonal(d::SpikedWishart)
    a = d.beta * d.n/2
    
    # diagonal and superdiagonal of B
    Bdv = [rand(Chi(2*a - d.beta*k))/sqrt(d.beta) for k in 0:(d.p-1)]
    Bev = [rand(Chi(d.beta*(d.p-k)))/sqrt(d.beta) for k in 1:(d.p-1)]

    if length(d.spikes) == 1
        Bdv[1] *= sqrt(1 + d.spikes[1])
    end
    
    # diagonal and off-diagonals of BB'
    dv = Bdv.^2 + [0; Bev.^2]
    ev = Bdv[1:end-1] .* Bev
    
    SymTridiagonal(dv, ev) * scaling(d)
end

# Multi-spike version only implemented for beta=1
function randbanded(d::SpikedWishart)
    @assert d.beta in [1, 2]
    r = length(d.spikes)

    if d.beta == 1
        U = BandedMatrix{Float64}(undef, (d.p, d.p), (0, r))
    elseif d.beta == 2
        U = BandedMatrix{Complex{Float64}}(undef, (d.p, d.p), (0, r))
    end

    dv = [rand(Chi(d.beta*(d.n - k + 1)))/sqrt(d.beta) for k in 1:d.p]
    @. dv[1:r] *= sqrt(1 + d.spikes)

    U[band(0)] .= dv
    
    for k = 1:(r-1)
        if d.beta == 1
            ev = randn(d.p - k)
        elseif d.beta == 2
            ev = (randn(d.p - k) + im * randn(d.p-k))/sqrt(2)
        end
        
        @. ev[1:(r-k)] *= sqrt(1 + d.spikes[(k+1):end])

        U[band(k)] .= ev
    end

    U[band(r)] .= [rand(Chi(d.beta*(d.p - k)))/sqrt(d.beta) for k in r:(d.p-1)]

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
    sigma = @. l * sqrt(gamma * (2/d.beta) * (1 - gamma/cspikes^2)) / sqrt(d.n)

    if d.scaled == false
        mu *= d.n
        sigma *= d.n
    end
    
    MvNormal(mu, Diagonal(sigma .^ 2))
end
