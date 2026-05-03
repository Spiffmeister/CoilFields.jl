# CoilFields.jl

A package for doing things with coils for fusion plasmas.


```
julia> using CoilFields
```

## Reading in coil sets

```julia
coilset = readcoilset("test/coilset")
```


There are inbuilt plotting recipes using [Makie]() for plotting `CoilSet` and single `Coil` objects
```julia
using GLMakie
plotcoils(coilset)
```

## Computing the magnetic field

The magnetic field can be evaluated by calling the `biot_savart` function,

```julia
pt = zeros(3)
B = biot_savart(coilset, pt, CompactLinear)
```

We can also evaluate at a number of points by handing a vector of points to the `biot_savart` function, note that we also export the in-place variants (called with the `!`).

```julia
pts = [rand(3) for _ in 1:100]
B = [zeros(3) for _ in eachindex(pts)]
biot_savart!(B, coilset, pts, CompactLinear)
```
