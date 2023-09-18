
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
    W=MathieuWron(ν,C_2k,index)
    return a,C_2k,index,W
end

"""
`W=MathieuWron(ν,q)` or `W=MathieuWron(ν,C_k,index)`.

Return the Wronskian.

For Mathieu's equation

y'' + (a- 2 q cos( 2z )) y = 0,

the Wronskian is defined by 

```math
\\frac{\\dot{f} f^*-f \\dot{f^*}}{2i}
```

where ``f=e^{i\\nu z}\\sum_k{C_{k}e^{i2kz}}`` with ``\\sum_k{C_{k}^2}=1`` is the solution of the Mathieu equation and ``f^*`` is its complex conjugate. 
"""
function MathieuWron(ν,q)
    _,C_2k,index=MathieuCharVecλ(ν,q)
    return MathieuWron(ν,C_2k,index)
end

function MathieuWron(ν,C_2k::Vector,index::Int)
    W=0.0
    for i in eachindex(C_2k),j in eachindex(C_2k)
        W+=C_2k[i]*C_2k[j]*(ν+(i-index+j-index))
    end
    return W # For use in real physical motion, W should be multiplied by ω_d/2 with normalization satisfying f(0)=1.
end

"""
``f=e^{i\\nu z}\\sum_k{C_{2k}e^{i2kz}}`` with ``\\sum_k{C_{2k}^2}=1`` is the solution of the Mathieu equation.
"""
function MathieuFunc(ν,q,z)
    _,C_2k,index=MathieuCharVecλ(ν,q)
    f=sum(C_2k.*exp.(im*(2*collect(1-index:length(C_2k)-index).+ν)*z))
    # TODO use @inbound for to speedup
    return f
end

"""
``\\partial f/\\partial z``
"""
function MathieuFuncPrime(ν,q,z)
    _,C_2k,index=MathieuCharVecλ(ν,q)
    f=sum(im*(2*collect(1-index:length(C_2k)-index).+ν).*C_2k.*exp.(im*(2*collect(1-index:length(C_2k)-index).+ν)*z))
    # TODO use @inbound for to speedup
    return f
end


""" 
[ν,c]=mathieu_mathieuexp(a,q;ndet::Int=20)

This program evaluates the characteristic exponent ν, 
corresponding to solutions of Mathieu Equation 
y''(t)+(a-2q cos(2t)) y(t)=0;
where a and q  are fixed real variables.

ndet is a positive integer number: it is the matrix dimension used 
in the algorithm. Precision increases by increasing ndet. 
Default value is ndet=20

The alghoritm consider two different cases: 
a=(2k)^2 or not (k integer).
ν is such that its real part belongs to the interval [0,2]
Of course, every other solutions are obtained by the formula
±ν+2k, with k integer.

"""
function MathieuExponent(a,q;ndet::Int=20,has_img::Bool=true,max_ndet::Int=1000)
    x=(a>=0)&& sqrt(abs(a))/2%1==0
    N=2*ndet+1 #matrix dimension
    a,q=float.(promote(a,q))
    d=q./((2*(-ndet:ndet) .+x).^2 .-a)
    m=Tridiagonal(d[2:N], ones(N), d[1:N-1])
    delta=det(m)
    if x
        if 0<=delta<=1
            alpha=acos(2*delta-1)/pi
            ν=mod(alpha,2) #modular reduction to the solution [0,2]
            H_nu=SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1))
            ck=eigvecs(H_nu,[0.0])[:,1]
            return ν,ck
        elseif has_img==true
            alpha=acos(2*Complex(delta)-1)/pi
            ν=alpha*(2*(imag(alpha)>=0)-1) #change an overall sign so that the imaginary part is always positive.
            ν=mod(real(ν),2)+im*imag(ν) #modular reduction to the solution [0,2]
            q=Complex(q)
            H_nu=Matrix(SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1)))
            vals,vecs=eigen(H_nu)
            _,idx=findmin(abs,vals)
            return ν,vecs[:,idx]
        elseif ndet<max_ndet
            MathieuExponent(a,q;ndet=2*ndet,has_img=false,max_ndet=max_ndet)
        else
            @warn "Expect real output for a=$a and q=$q, but the result is complex even for ndet=$ndet."
            alpha=acos(2*Complex(delta)-1)/pi
            ν=alpha*(2*(imag(alpha)>=0)-1) #change an overall sign so that the imaginary part is always positive.
            ν=mod(real(ν),2)+im*imag(ν) #modular reduction to the solution [0,2]
            q=Complex(q)
            H_nu=Matrix(SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1)))
            vals,vecs=eigen(H_nu)
            _,idx=findmin(abs,vals)
            return ν,vecs[:,idx]
        end
    else
        beta=delta*sin(pi*sqrt(Complex(a))/2)^2
        beta=real(beta) # beta should be real.
        if 0<=beta<=1
            alpha=2*asin(sqrt(beta))/pi
            ν=mod(alpha,2) #modular reduction to the solution [0,2]
            H_nu=SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1))
            ck=eigvecs(H_nu,[0.0])[:,1]
            return ν,ck
        elseif has_img==true
            alpha=2*asin(sqrt(Complex(beta)))/pi
            ν=alpha*(2*(imag(alpha)>=0)-1) #change an overall sign so that the imaginary part is always positive.
            ν=mod(real(ν),2)+im*imag(ν) #modular reduction to the solution [0,2]
            q=Complex(q)
            H_nu=Matrix(SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1)))
            vals,vecs=eigen(H_nu)
            _,idx=findmin(abs,vals)
            return ν,vecs[:,idx]
        elseif ndet<max_ndet
            MathieuExponent(a,q;ndet=min(2*ndet,max_ndet),has_img=false,max_ndet=max_ndet)
        else
            @warn "Expect real output for a=$a and q=$q, but the result is complex even for ndet=$ndet."
            alpha=2*asin(sqrt(Complex(beta)))/pi
            ν=alpha*(2*(imag(alpha)>=0)-1) #change an overall sign so that the imaginary part is always positive.
            ν=mod(real(ν),2)+im*imag(ν) #modular reduction to the solution [0,2]
            q=Complex(q)
            H_nu=Matrix(SymTridiagonal((ν .+ 2*(-ndet:ndet)).^2 .- a,q*ones(N-1)))
            vals,vecs=eigen(H_nu)
            _,idx=findmin(abs,vals)
            return ν,vecs[:,idx]
        end
    end
end
