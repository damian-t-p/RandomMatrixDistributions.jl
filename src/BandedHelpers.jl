export eigmax

function eigmax(A::Symmetric{T,<:BandedMatrix{T}}) where T <: Real
    maximum(KrylovKit.eigsolve(A, issymmetric=true)[1])
end

function eigmax(A::Hermitian{T,<:BandedMatrix{T}}) where T
    maximum(KrylovKit.eigsolve(A, ishermitian=true)[1])
end
