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
        throw(">1 spike not yet implemented for Wigner matrices")
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



