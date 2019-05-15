export Wachter

"""
Wachter Distribution

Let Σ<sub>1</sub> be Σ<sub>2</sub> pxp covariance matrices of n<sub>1</sub> and n<sub>2</sub> standard normal observations respectively.

If p/n<sub>1</sub> → γ<sub>1</sub> and p/n<sub>2</sub> → γ<sub>2</sub>, then Σ<sub>1</sub>Σ<sub>2</sub><sup>-1</sup> has a limiting spectral
distribution of Wachter(γ<sub>1</sub>, γ<sub>2</sub>).
"""
struct Wachter <: ContinuousUnivariateDistribution
    gamma1::Real
    gamma2::Real
    Wachter(gamma1, gamma2) = 0 < gamma2 < 1 ? new(gamma1, gamma2) : error("Gamma2 must be in (0, 1)")
end

function pdf(d::Wachter, x)
    a = ((1 - sqrt(d.gamma1 + d.gamma2 - d.gamma1 * d.gamma2))/(1 - d.gamma2))^2
    b = ((1 + sqrt(d.gamma1 + d.gamma2 - d.gamma1 * d.gamma2))/(1 - d.gamma2))^2

    if a < x < b
        return (1 - d.gamma2) * sqrt((b - x) * (x - a))/(2 * pi * x * (d.gamma1 + d.gamma2 * x))
    else
        return 0
    end
end
