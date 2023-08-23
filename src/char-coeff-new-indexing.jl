
"""
MathieuCharA(ν,q)

char value A_ν for Mathieu's equation

y'' + (A_ν - 2 q cos( 2z )) y = 0

where

q ∈ ℝ - parameter
ν ∈ ℝ - characteristic exponent

"""
function MathieuCharA(ν::Real,q::Real)
    if isinteger(ν)
        return MathieuCharA(Int(ν),q)
    else
        return MathieuCharλ(ν,q)
    end
end

function MathieuCharA(ν::Int,q::Real)
    if ν==0
        return charλ(abs(q), 0.0; k=1:1)[1]
    else
        return charA(q; k=abs(ν):abs(ν))[1]
    end
end

"""
MathieuCharB(ν,q)

char value B_ν for Mathieu's equation

y'' + (B_ν - 2 q cos( 2z )) y = 0

where

q ∈ ℝ       - parameter
ν ∈ ℝ  - fractional part of the non-integer order

"""
function MathieuCharB(ν::Real,q::Real)
    if isinteger(ν)
        return MathieuCharB(Int(ν),q)
    else
        return MathieuCharλ(ν,q)
    end
end

function MathieuCharB(ν::Int,q::Real)
    return charB(q; k=abs(ν):abs(ν))[1]
end


function MathieuCharλ(ν::Real,q::Real) # reduced = true
    #nu = reduced ? rem(nu_+1,2)-1 : nu_;
    ν_abs=abs(ν)
    (k,ν_)=divrem(ν_abs,2,RoundNearest)

    # Set matrix size using formula from Shirts paper (1993), Eqs. (2.1)-(2.2).
    # nu0 = nu + maximum(k)
    C = (8.46 + 0.444*ν)/(1 + 0.085*ν)
    D =  (0.24 + 0.0214*ν)/(1 + 0.059*ν)
    N = ceil(Int, (ν + 2 + C*abs(q)^D)/2) # matrix size is 2N+1

    (two, q2, nu2) = float.(promote(2, q, ν_)) # d0 and d1 must be of the same type
    d0 = (two .* (-N:N) .+ nu2).^2
    d1 = fill(q2, 2 * N)
    A = SymTridiagonal(d0, d1)
    a = eigvals!(A, trunc(Int,ν)+1:trunc(Int,ν)+1)[1]
    return a
end

"""
Return eigenvalue, eigenvector (Fourier coefficient) of the Mathieu characteristic value problem and the index of the k=0 Fourier coefficient.
"""
function MathieuCharVecλ(ν::Real,q::Real) # reduced = true
    #nu = reduced ? rem(nu_+1,2)-1 : nu_;
    ν_abs=abs(ν)
    (k,ν_)=divrem(ν_abs,2,RoundNearest)

    # Set matrix size using formula from Shirts paper (1993), Eqs. (2.1)-(2.2).
    # nu0 = nu + maximum(k)
    C = (8.46 + 0.444*ν)/(1 + 0.085*ν)
    D =  (0.24 + 0.0214*ν)/(1 + 0.059*ν)
    N = ceil(Int, (ν + 2 + C*abs(q)^D)/2) # matrix size is 2N+1

    (two, q2, nu2) = float.(promote(2, q, ν_)) # d0 and d1 must be of the same type
    d0 = (two .* (-N:N) .+ nu2).^2
    d1 = fill(q2, 2 * N)
    A = SymTridiagonal(d0, d1)
    vals,vecs = eigen!(A)
    center_index=Int(N+1+k)
    return vals[trunc(Int,ν)+1],vecs[:,trunc(Int,ν)+1],center_index
end

"""
Return eigenvalue, eigenvector (Fourier coefficient) of the Mathieu characteristic value problem, the index of the k=0 Fourier coefficient and the Wronskian.

For Mathieu's equation

y'' + (B_ν - 2 q cos( 2z )) y = 0,

the Wronskian is defined by 

```math
\\frac{\\dot{f} f^*-f \\dot{f^*}}{2i}
```

where ``f`` is the solution of the Mathieu equation and ``f^*`` is its complex conjugate. 
"""
function MathieuCharVecWronλ(ν::Real,q::Real)
    a,C_2k,index=MathieuCharVecλ(ν,q)
    W=0.0
    for i in eachindex(C_2k),j in eachindex(C_2k)
        W+=C_2k[i]*C_2k[j]*(ν+(i-index+j-index))
    end
    return a,C_2k,index,W
end

"""
Return the Wronskian.

For Mathieu's equation

y'' + (B_ν - 2 q cos( 2z )) y = 0,

the Wronskian is defined by 

```math
\\frac{\\dot{f} f^*-f \\dot{f^*}}{2i}
```

where ``f`` is the solution of the Mathieu equation and ``f^*`` is its complex conjugate. 
"""
function MathieuWron(ν,q)
    _,C_2k,index=MathieuCharVecλ(ν,q)
    W=0.0
    for i in eachindex(C_2k),j in eachindex(C_2k)
        W+=C_2k[i]*C_2k[j]*(ν+(i-index+j-index))
    end
    return W
end


