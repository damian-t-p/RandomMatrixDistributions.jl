export Jacobi

"""
Jacobi Distribution

If E ~ Wishart<sub>p</sub>(I, n<sub>1</sub>) and H ~ Wishart<sub>p</sub>(I, n<sub>2</sub>) are independent,
then E(E + H)<sup>-1</sup> has Jacobi(n<sub>1</sub>, n<sub>2</sub>, p) distribution.

If λ<sub>i</sub> are the eigenvalues of EH<sup>-1</sup> and μ<sub>i</sub> are the Jacobi eigenvalues,
then μ<sub>i</sub> = λ<sub>i</sub>/(1 + λ<sub>i</sub>) and  λ<sub>i</sub> = μ<sub>i</sub>/(1 - μ<sub>i</sub>)
"""
struct Jacobi <: ContinuousMatrixDistribution
    beta::Integer
    n1::Integer
    n2::Integer
    p::Integer
end

function randreduced(d::Jacobi)
    randtridiagonal(d)
end

function randtridiagonal(d::Jacobi)
    
    a = d.beta/2 * (d.n1 - d.p + 1) - 1;
    b = d.beta/2 * (d.n2 - d.p + 1) - 1;
    
    alphaeven = [2*rand(Beta((2*d.p - k - 2)*d.beta/4 + a + 1, (2*d.p - k - 2) * d.beta/4 + b + 1)) - 1 for k in 0:(2 * d.p-2) if k%2 == 0]
    alphaodd = [2*rand(Beta((2*d.p - k - 1)*d.beta/4, (2*d.p - k - 3) * d.beta/4 + a + b + 2)) - 1 for k in 0:(2 * d.p-2) if k%2 == 1]

    alphaevenleft = [alphaeven; 0]
    alphaevenright = [0; alphaeven]
    alphaodd = [-1; alphaodd; -1]  
    
    dv = (1 .- alphaodd) .* alphaevenleft - (1 .+ alphaodd).* alphaevenright
    dv = dv[1:(end-1)]

    ev = sqrt.((1 .- alphaodd[1:(end-1)]) .* (1 .- alphaeven.^2) .* (1 .+ alphaodd[2:end]))
    ev = ev[1:end-1]
                            
    (SymTridiagonal(dv, ev) + 2I)/4
end
