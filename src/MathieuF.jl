module MathieuF
    using LinearAlgebra
    using MathieuFunctions

    export  MathieuCharA,MathieuCharB

    include("char-new-indexing.jl")
    include("fourier-coeffs.jl")
end