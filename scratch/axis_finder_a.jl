using CoilFields

using NonlinearSolve
using OrdinaryDiffEq
using LinearAlgebra

using GLMakie

coilset = ReadCoilSet("./test/coilset", skipstart=3)





"""
Our intersection plane Σ is a circle centered at C_0 with a normal vector n_c and radius r_c
The circle is assumed to be large enough such that any point which passes though it will
return to it in finite time. I.e. the DE has the property that,
f : Σ → Σ

"""
function event(u, t, integrator)
    # C₀ = [10.0, 0.0, 0.0]
    # circle_normal = [1.0, 0.0, 0.0]
    # circle_radius = 1.0

    # uᵢ = integrator.uprev
    # First check if the point crosses the plane defined by the normal vector
    # dot(u - circle_normal, uᵢ - circle_normal) # < 0 if crossed
    # (dot(u, circle_normal) - circle_radius) * circle_normal
    # dot(uᵢ, circle_normal)

    # Define a normal vector to a plane, event triggered when h(u_i) * h(u_{i-1}) < 0 where h(u) is the plane
    normal = [1.0, 0.0, 0.0]
    ret = one(eltype(u))
    # If the plane has been crossed, then if sign(u) < 0 then sign(u) h(u_i) h(u_{i-1}) > 0
    # else sign(u) h(u_i) h(u_{i-1}) < 0 so event triggered
    ret = dot(u, normal)
    # sign((u - integrator.uprev, normal))
    # ret = sign(dot(normal, u)) * dot(u, normal)
    # end
    return ret
end


event_affect!(integrator) = terminate!(integrator)
event_affect_neg! = nothing

cb = ContinuousCallback(event, event_affect!, (affect_neg!)=event_affect_neg!, save_positions=(true, false))

x₀ = [10.0, 0.0, 0.0]
fieldline_ode = ODEProblem((ẋ, x, p, t) -> CoilFields.field_line!(ẋ, x, p, t, coilset), x₀, (0.0, 800))
# fieldline = solve(fieldline_ode, Tsit5(), callback=cb, save_everystep=false)
fieldline = solve(fieldline_ode, Tsit5())


"""
Find the point ``X`` such that ``X - x=0`` where ``x`` is determined by following a field line for one turn
"""
function axis_diff(X, coilset)
    fieldline = ODEProblem((ẋ, x, p, t) -> CoilFields.field_line!(ẋ, t, x, p, coilset), X, (0.0, 800))
    x = solve(fieldline).u
    x = X .- x
end


function find_axis(X, coilset)
    nlp = NonlinearProblem((X, p) -> axis_diff(X, coilset), X)
    nls = solve(nlp)
    return nls.u
end


us = hcat(fieldline.u...)

scatter(us[1, :], us[2, :], us[3, :])
