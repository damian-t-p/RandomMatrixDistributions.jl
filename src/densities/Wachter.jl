export Wachter

"""
    Wachter(γ₁, γ₂)

Wachter Distribution, where 0 ≤ `γ₂` < 1.

Let Σ₁ and Σ₂ be p×p covariance matrices of n₁ and n₂ standard normal observations respectively.

If p/n₁ → γ₁ and p/n₂ → γ₂, then Σ₁ Σ₂⁻¹ has a limiting spectral
distribution of `Wachter(γ₁, γ₂)`.
"""
struct Wachter <: ContinuousUnivariateDistribution
    gamma1::Real
    gamma2::Real
    Wachter(gamma1, gamma2) = 0 <= gamma2 < 1 ? new(gamma1, gamma2) : error("Gamma2 must be in [0, 1)")
end

function minimum(d::Wachter)
    ((1 - sqrt(d.gamma1 + d.gamma2 - d.gamma1 * d.gamma2))/(1 - d.gamma2))^2
end

function maximum(d::Wachter)
    ((1 + sqrt(d.gamma1 + d.gamma2 - d.gamma1 * d.gamma2))/(1 - d.gamma2))^2
end

function pdf(d::Wachter, x::Real)
    a = minimum(d)
    b = maximum(d)

    if a < x < b
        return (1 - d.gamma2) * sqrt((b - x) * (x - a))/(2 * pi * x * (d.gamma1 + d.gamma2 * x))
    else
        return 0
    end
end
