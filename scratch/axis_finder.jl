using CoilFields

using NonlinearSolve
using OrdinaryDiffEq
using LinearAlgebra

coilset = ReadCoilSet("./test/coilset", skipstart=3)





"""
Our intersection plane Σ is a circle centered at C_0 with a normal vector n_c and radius r_c
The circle is assumed to be large enough such that any point which passes though it will
return to it in finite time. I.e. the DE has the property that,
f : Σ → Σ

h(x) = {s = F(x,y,z) | σ(x,y,z) = 0}
"""
function event(u, t, integrator)
    # C₀ = [10.0, 0.0, 0.0]
    plane_normal = [0.0, 1.0, 0.0]
    # circle_radius = 1.0

    # uᵢ = integrator.uprev
    # First check if the point crosses the plane defined by the normal vector
    # dot(u - circle_normal, uᵢ - circle_normal) # < 0 if crossed
    # Project the outgoing point onto the plane, then check if this is within the radius
    # proj_u = u - dot(u, circle_normal) / norm(circle_normal) * circle_normal
    # norm(proj_u - C₀) < circle_radius # True if within radius
    # Need to turn this into a condition
    dot(u, plane_normal)
end


event_affect!(integrator) = terminate!(integrator)
affect_neg!(integrator) = nothing

cb = ContinuousCallback(event, event_affect!; (affect_neg!)=nothing, save_positions=(true, false))

"""
Find the point ``X`` such that ``X - x=0`` where ``x`` is determined by following a field line for one turn
"""
function axis_diff(X, coilset)
    fieldline = ODEProblem((ẋ, x, p, t) -> CoilFields.field_line!(ẋ, x, p, t, coilset), [X[1], 0.0, X[2]], (0.0, 800))
    sol = solve(fieldline, callback=cb, save_everystep=false)
    x = sol.u[end]
    return [X[1] - x[1], X[2] - x[3]]
end


function find_axis(X, coilset)
    nlp = NonlinearProblem((X, p) -> axis_diff(X, coilset), X)
    nls = solve(nlp)
    return [nls.u[1], 0.0, nls.u[2]]
end



X = [10.0, 0.0]
# x = axis_diff(X, coilset)

x_ax = find_axis(X, coilset)



fieldline = ODEProblem((ẋ, x, p, t) -> CoilFields.field_line!(ẋ, x, p, t, coilset), [10.0, 0.0, 0.0], (0.0, 800))
sol = solve(fieldline, Tsit5(), callback=cb, save_everystep=false)




# using GLMakie
# u = sol.u...
# scatter(hcat(sol.u...))
