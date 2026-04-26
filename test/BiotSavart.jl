using CoilFields
using Test

using SpecialFunctions: ellipk, ellipe
using LinearAlgebra
using StaticArrays


μ_0 = CoilFields.μ₀

# Define a single magnetic coil







# A = CoilFields.biot_savart_A(circular_coil, pts, CoilFields.CompactLinear())






# Define the circular coil
npts = 10_000
J = 100.0
coil_points = Tuple(SVector(cos(θ), sin(θ), 0.0) for θ in range(0.0, 2π, npts))
circular_coil = Coil(coil_points, J, npts)

# Define sampling points
pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
origin = zeros(3)

circular_coil_axis_mod_B(z, R, J) = μ_0 * J * R^2 / (2 * (z^2 + R^2)^(3 / 2))


@testset "Compact linear segment testing" begin
    # Single evaluation of the Boit-Savart at the coil origin

    @testset "Single point on axis" begin
        B = biot_savart(circular_coil, origin, CompactLinear())
        B_exact = circular_coil_axis_mod_B(zero(eltype(origin)), one(eltype(origin)), J)
        @test (norm(B) - B_exact) ≤ 1e-10
    end

    @testset "Multiple points along axis" begin
        B = map(x -> biot_savart(circular_coil, x, CompactLinear()), pts)
    end

end






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


# Single point on coil
biot_savart(circular_coil, [0.0, 0.0, 0.0], CompactLinear())
B_circ_exact = circular_coil_axis_mod_B(0.0, one(eltype(pts[1])), J)


coilset = CoilSet([circular_coil, circular_coil])
biot_savart(coilset, [0.0, 0.0, 0.0], CompactLinear())
2 * B_circ_exact


B_circ_exact = circular_coil_axis_mod_B.(last.(pts), one(eltype(pts[1])), J)
