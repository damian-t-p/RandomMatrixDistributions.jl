export SpikedWishart

struct SpikedWishart <: ContinuousMatrixDistribution
    beta::Integer
    n::Integer
    p::Integer
    spikes::Array{Float64, 1}
end

SpikedWishart(beta, n, p) = SpikedWishart(beta, n, p, [])

# SAMPLERS

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
    Bdv = [rand(Chi(2*a - d.beta*k)) for k in 0:(d.p-1)]
    Bev = [rand(Chi(d.beta*(d.p-k))) for k in 1:(d.p-1)]

    if length(d.spikes) == 1
        Bdv[1] *= sqrt(1 + d.spikes[1])
    end
    
    # diagonal and off-diagonals of BB'
    dv = Bdv.^2 + [0; Bev.^2]
    ev = Bdv[1:end-1] .* Bev

    SymTridiagonal(dv, ev)
end

# Multi-spike version only implemented for beta=1
function randbanded(d::SpikedWishart)
    @assert d.beta == 1
    r = length(d.spikes)

    U = BandedMatrix{Float64}(undef, (d.p, d.p), (0, r))

    dv = [rand(Chi(d.n - k + 1)) for k in 1:d.p]
    @. dv[1:r] *= sqrt(1 + d.spikes)

    U[band(0)] .= dv
    
    for k = 1:(r-1)
        ev = randn(d.p - k)
        @. ev[1:(r-k)] *= sqrt(1 + d.spikes[(k+1):end])

        U[band(k)] .= ev
    end

    U[band(r)] .= [rand(Chi(d.p - k + 1)) for k in 1:(d.p-r)]

    Symmetric(U' * U)
end
