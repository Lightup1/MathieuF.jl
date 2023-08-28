# MathieuF.jl

Julia package for Mathieu Functions with function forms similar to Mathieu related functions in Mathematica.

Mathieu functions are the eigenvalues and eigenfunction of *Mathieu's
equation* (for more details, see [NIST's Digital Library of
Mathematical Functions](http://dlmf.nist.gov/28)).

Related Package: [MathieuFunctions.jl](https://github.com/BBN-Q/MathieuFunctions.jl)

## Highlights
- Supporting output of Fourier coefficients of non-integer order eigenfunctions
- Supporting output of the related Wronskian.

## TODO
- support a,q as input using continued fraction, see [the algorithm](https://www.jstor.org/stable/2003814).
- Adopting MTIEU2 algorithm for faster speed when only single eigenvalue and the corresponding eigenfunction are needed, see [Shirts](http://dl.acm.org/citation.cfm?id=155796).
- Accuracy problem for large q or order.
