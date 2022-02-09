export Jacobi

"""
    Jacobi(β::Int, n₁::Int, n₂::Int, p::Int)

Distribution of a p×p Jacobi matrix.

If E ~ Wishartₚ(I, n₁) and H ~ Wishartₚ(I, n₂) are independent with Dyson parameter β,
then E(E + H)⁻¹ has `Jacobi(β, n₁, n₂, p)` distribution.

If λᵢ are the eigenvalues of EH⁻¹ and μᵢ are the Jacobi eigenvalues,
then μᵢ = λᵢ/(1 + λᵢ) and λᵢ = μᵢ/(1 - μᵢ).
"""
struct Jacobi <: ContinuousMatrixDistribution
    beta::Integer
    n1::Integer
    n2::Integer
    p::Integer
end

# PROPERTIES

Base.size(d::Jacobi) = (d.p, d.p)

# SAMPLERS

function randreduced(rng::AbstractRNG, d::Jacobi)
    randtridiagonal(rng, d)
end

function randtridiagonal(rng::AbstractRNG, d::Jacobi)
    
    a = d.beta/2 * (d.n1 - d.p + 1) - 1;
    b = d.beta/2 * (d.n2 - d.p + 1) - 1;
    
    alphaeven = [2*rand(rng, Beta((2*d.p - k - 2)*d.beta/4 + a + 1, (2*d.p - k - 2) * d.beta/4 + b + 1)) - 1 for k in 0:(2 * d.p-2) if k%2 == 0]
    alphaodd = [2*rand(rng, Beta((2*d.p - k - 1)*d.beta/4, (2*d.p - k - 3) * d.beta/4 + a + b + 2)) - 1 for k in 0:(2 * d.p-2) if k%2 == 1]

    alphaevenleft = [alphaeven; 0]
    alphaevenright = [0; alphaeven]
    alphaodd = [-1; alphaodd; -1]  
    
    dv = (1 .- alphaodd) .* alphaevenleft - (1 .+ alphaodd).* alphaevenright
    dv = dv[1:(end-1)]

    ev = sqrt.((1 .- alphaodd[1:(end-1)]) .* (1 .- alphaeven.^2) .* (1 .+ alphaodd[2:end]))
    ev = ev[1:end-1]
                            
    (SymTridiagonal(dv, ev) + 2I)/4
end
