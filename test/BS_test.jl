using Revise
using CoilFields
using StaticArrays

# using ProfileView, Cthulhu
using BenchmarkTools

npts = 10_000
J = 100.0
coil_points = Tuple(SVector(cos(θ), sin(θ), 0.0) for θ in range(0.0, 2π, npts))
circular_coil = CoilFields.Coil(coil_points, J, npts)


# @benchmark CoilFields.biot_savart($circular_coil, $[0.0, 0.0, 0.0], $CoilFields.CompactLinear())

pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
B = [zeros(3) for _ in eachindex(pts)];
@benchmark CoilFields.biot_savart!($B, $circular_coil, $pts, $CoilFields.CompactLinear())
