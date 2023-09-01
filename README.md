# MathieuF.jl

Julia package for Mathieu Functions with function forms similar to Mathieu related functions in Mathematica.

Mathieu functions are the eigenvalues and eigenfunction of *Mathieu's
equation* (for more details, see [NIST's Digital Library of
Mathematical Functions](http://dlmf.nist.gov/28)).

Related Package: [MathieuFunctions.jl](https://github.com/BBN-Q/MathieuFunctions.jl)

## Highlights
- Support Fourier coefficients of non-integer order eigenfunctions
- Support a,q as input, see [Mathieu functions toolbox](https://atoms.scilab.org/toolboxes/Mathieu/4.0.61)
- Support output of the related Wronskian.


## Examples

```julia
nu,ck=MathieuExponent(a,q)
```
where `nu` is the characteristic exponent and vector `ck` is the Fourier coefficients of the eigenfunction with `norm(ck)≈1`.
Note that `nu` is reduced to the interval `[0,2]` and `c0` corresponds to `ck[(length(ck)-1)÷2]` with the reduced `nu`.

```julia
W=MathieuWron(nu,ck::Vector,index::Int)
```
where `W` is the Wronskian of the eigenfunction with and `index` refers to the index of `c0` in `ck`. 
For example, 
```julia
a=0.1;q=0.5;
nu,ck=MathieuExponent(a,q)
idx=(length(ck)-1)÷2+1
W=MathieuWron(nu,ck,idx)
```

If one knows `nu` (not reduced) and `q`, one can use
```julia
W=MathieuWron(nu,q)
```


