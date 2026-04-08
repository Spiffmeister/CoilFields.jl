# using Revise
using CoilFields

using Profile

using SpecialFunctions: ellipk, ellipe
using LinearAlgebra
using StaticArrays


import PhysicalConstants.CODATA2022: μ_0


# Define a single magnetic coil
npts = 10_000

coil_points = Tuple(SVector(cos(θ), sin(θ), 0.0) for θ in range(0.0, 2π, npts))

J = 100.0
circular_coil = CoilFields.Coil(coil_points, J, npts)


# Single evaluation of the Boit-Savart at the coil origin
CoilFields.Biot_Savart(circular_coil, [0.0, 0.0, 0.0], CoilFields.CompactLinear())


# Define points along the z-axis
pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
B = [zeros(3) for _ in eachindex(pts)];

B = CoilFields.Biot_Savart(circular_coil, pts, CoilFields.CompactLinear())
CoilFields.Biot_Savart!(B, circular_coil, pts, CoilFields.CompactLinear())


A = CoilFields.Biot_Savart_A(circular_coil, pts, CoilFields.CompactLinear())


# using BenchmarkTools

# @benchmark CoilFields.Biot_Savart!($zeros(3),$circular_coil, $[0.0,0.0,0.0], $CoilFields.CompactLinear)
# @benchmark CoilFields.Biot_Savart($circular_coil_s, $[0.0,0.0,0.0], $CoilFields.CompactLinear)
# @benchmark CoilFields.Biot_Savart!($zeros(3),$circular_coil_s, $[0.0,0.0,0.0], $CoilFields.CompactLinear)

# @benchmark CoilFields.Biot_Savart($circular_coil_s,$pts,$CoilFields.CompactLinear)
# @benchmark CoilFields.Biot_Savart!($B,$circular_coil_s,$pts,$CoilFields.CompactLinear)


# using Profile

# CoilFields.Biot_Savart!(B,circular_coil,pts,CoilFields.CompactLinear)
# Profile.clear_malloc_data()
# CoilFields.Biot_Savart!(B,circular_coil,pts,CoilFields.CompactLinear)


# using ProfileView
# using Cthulhu
# @profview CoilFields.Biot_Savart!(B,circular_coil_s,pts,CoilFields.CompactLinear)
# @profview CoilFields.Biot_Savart!(B,circular_coil_s,pts,CoilFields.CompactLinear)


circular_coil_axis_mod_B(z, R, J) = μ_0 * J * R^2 / (2 * (z^2 + R^2)^(3 / 2))




function analytic_vector_potential()
    val1 = z^2 + (R + ρ)^2
    val2 = √(4 * R * ρ / val1)

    A = √(val1) / (2 * R * ρ) * (
        (1 - (2 * R * ρ) / val1) * ellipk(val2) - ellipe(val2)
    )

    A .= I * μ_0 / 2π * A

    A_cartesian = norm(A) * []

    return A
end




B_circ_exact = circular_coil_axis_mod_B.(last.(pts), one(eltype(pts[1])), J)
