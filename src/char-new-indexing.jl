export  MathieuCharA, 
        MathieuCharB


"""
MathieuCharA(ν,q)

char value A_ν for Mathieu's equation

y'' + (A_ν - 2 q cos( 2z )) y = 0

where

q ∈ ℝ       - parameter
ν ∈ ℝ  - fractional part of the non-integer order

"""
function MathieuCharA(ν::Real,q::Real)
    if isinteger(ν)
        return MathieuCharA(Int(ν),q)
    else
        v_trunc=trunc(ν)
        return charλ(q, ν-v_trunc; k=v_trunc:v_trunc)
    end
end

function MathieuCharA(ν::Int,q::Real)
    charA(q; k=ν:ν)
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
        return charλ(q, ν-v_trunc; k=v_trunc:v_trunc)
    end
end

function MathieuCharB(ν::Int,q::Real)
    charB(q; k=ν:ν)
end