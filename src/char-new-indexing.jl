
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
        ν_abs=abs(ν)
        ν_abs_trunc=trunc(Int,ν_abs)
        return charλ(q, ν_abs-ν_abs_trunc; k=ν_abs_trunc+1:ν_abs_trunc+1)[1]
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
        ν_abs=abs(ν)
        ν_abs_trunc=trunc(Int,ν_abs)
        return charλ(q, ν_abs-ν_abs_trunc; k=ν_abs_trunc+1:ν_abs_trunc+1)[1]
    end
end

function MathieuCharB(ν::Int,q::Real)
    return charB(q; k=abs(ν):abs(ν))[1]
end