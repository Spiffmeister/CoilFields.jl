using CoilFields

using NonlinearSolve
using OrdinaryDiffEq
using LinearAlgebra

coilset = ReadCoilSet("../test/coilset")





"""
Our intersection plane Σ is a circle centered at C_0 with a normal vector n_c and radius r_c
The circle is assumed to be large enough such that any point which passes though it will
return to it in finite time. I.e. the DE has the property that,
f : Σ → Σ

"""
function event(u, t, integrator, ζ)
    C₀ = [10.0, 0.0, 0.0]
    circle_normal = [1.0, 0.0, 0.0]
    circle_radius = 1.0

    uᵢ = integrator.uprev
    # First check if the point crosses the plane defined by the normal vector
    dot(u - circle_normal, uᵢ - circle_normal) # < 0 if crossed

    (dot(u, circle_normal) - circle_radius) * circle_normal
    dot(uᵢ, circle_normal)
end


event_affect!(integrator) = terminate!(integrator)

cb = ContinuousCallback(event, event_affect!, save_positions=(true, false))

"""
Find the point ``X`` such that ``X - x=0`` where ``x`` is determined by following a field line for one turn
"""
function axis_diff(X, coilset)
    fieldline = ODEProblem((ẋ, x, p, t) -> CoilFields.fieldline!(ẋ, t, x, p, coilset), X, (0.0, 800))
    x = solve(fieldline).u
    x = X .- x
end


function find_axis(X, coilset)
    nlp = NonlinearProblem((X, p) -> axis_diff(X, coilset), X)
    nls = solve(nlp)
    return nls.u
end


coilset()
