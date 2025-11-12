# CoilFields.jl

A package for doing things with coils for fusion plasmas.


```
julia> using CoilSet
```

## Reading in coil sets

```julia
coilset = ReadCoilSet("test/coilset")
```


There are inbuilt plotting recipes using [Makie]() for plotting `CoilSet` and single `Coil` objects
```julia
using GLMakie
plotcoils(coilset)
```

## Computing the magnetic field

The magnetic field can be evaluated by calling the `Biot_Savart` function,

```julia
pt = zeros(3)
B = Biot_Savart(coilset, pt, CompactLinear)
```

We can also evaluate at a number of points by handing a vector of points to the `Biot_Savart` function, note that we also export the in-place variants (called with the `!`).

```julia
pts = [rand(3) for _ in 1:100]
B = [zeros(3) for _ in eachindex(pts)]
Biot_Savart!(B, coilset, pts, CompactLinear)
```
