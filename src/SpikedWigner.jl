export SpikedWigner

struct SpikedWigner <: ContinuousMatrixDistribution
    beta::Integer
    n::Integer
    spikes::Array{Float64, 1}
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

