using CoilFields




coil_points = [[cos(θ), sin(θ), 0.0] for θ in range(0.0, 2π, npts)]
J = 100.0
circular_coil = CoilFields.Coil(coil_points, J, npts)



# In (r,θ,z) then on axis
# B = B_z = μ₀ I R² / 2 (z² + R²)^(3/2)
# Since B = ∇⨯A
# B_z = ∂_r A - ∂_θ A ⟹ A ≠ 0 on axis
# Since everything is symmetric about the point, ∂_θ A = 0


function analytic_vector_potential()
    val1 = z^2 + (R + ρ)^2
    val2 = √(4 * R * ρ / val1)

    √(val1) / (2 * R * ρ) * (
        (1 - (2 * R * ρ) / val1)
    )
end
