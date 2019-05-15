export eigmax

function eigmax(A::Symmetric{T,<:BandedMatrix{T}}) where T <: Real
    maximum(KrylovKit.eigsolve(A, issymmetric=true)[1])
end
