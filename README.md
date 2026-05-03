# CoilFields.jl

[![Documentation](https://shields.io/badge/docs-dev-blue)](https://spiffmeister.github.io/FaADE.jl/dev/)
[![Documentation](https://shields.io/badge/docs-dev-blue)](https://spiffmeister.github.io/CoilFields.jl/dev/)

A package for doing things with coils for fusion plasmas.

```
julia> using CoilFields
```

## Reading in coil sets

```julia
coilset = readcoilset("test/coilset")
```


There are inbuilt plotting recipes using [Makie]() for plotting `CoilSet` and single `Coil` objects,
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
