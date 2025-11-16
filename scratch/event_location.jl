using Revise

using GLMakie

using CoilFields

using OrdinaryDiffEq
using LinearAlgebra

coilset = ReadCoilSet("./test/coilset", skipstart=3)


function CartesianFieldLine!(ẋ, x, params, t)
    B = zeros(eltype(x), 3)
    Biot_Savart!(B, coilset, x, CompactLinear)

    ẋ .= B / norm(B)
end

prob = ODEProblem(CartesianFieldLine!, [10.0, 0.0, 0.0], (0.0, 200.0))
sol = solve(prob, Tsit5())

fl = FieldLine(sol.u, sol.t)


fcoil = plotcoils(coilset)

lines!(fcoil.axis, hcat(fl.position...)[1, :], hcat(fl.position...)[2, :], hcat(fl.position...)[3, :])

fcoil



# Consider the event defined by a thin ellipse
function event_ellipse(x, t, integrator, params)
    x₀, y₀, z₀, a, b, c = params
    h = (x[1] - x₀)^2 / a^2 + (x[2]^2 - y₀) / b^2 + (x[3] - z₀) / c^2 - one(eltype(x))
end
function event_ellipse_prime(x, params)
    x₀, y₀, z₀, a, b, c = params
    h′ = [
        2(x[1] - x₀) / a^2,
        2(x[2] - y₀) / b^2,
        2(x[3] - z₀) / c^2
    ]
end

# We then modify the field line ODE
# Redefine the normal field line ODE here and remove the in place for simplicity
function field_line_RHS!(x, p, t, coilset)
    B = zeros(eltype(x), 3)
    Biot_Savart!(B, coilset, x, CompactLinear)
    return B / norm(B)
end
function field_line_RHS_modified(x, p, t, coilset)
    field_line_RHS(x, p, t, coilset) / dot(event_ellipse_prime(x, p), field_line_RHS(x, p, t, coilset))
end

p = (10.0,)
