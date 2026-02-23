
struct FieldLine{TT}
    position::Vector{Vector{TT}}
    t::Vector{TT}
end

struct PoincarePlane{TT}
    points::Vector{Vector{TT}}
    ζ::TT
end

function field_line!(ẋ, x, p, t, coilset)
    B = zeros(eltype(x), 3)
    Biot_Savart!(B, coilset, x, CompactLinear)
    ẋ .= B / norm(B)
end


function field_line_branch!(ẋ, x, p, t, coilset)
    B = zeros(eltype(x), 3)
    Biot_Savart!(B, coilset, x, CompactLinear)
    for i in 1:3
        ẋ[i] = B[i] / norm(B)
    end
    ẋ[4] = atan(u[2], u[1]) - ζ
end



function _initialise_fieldlines(X₀::AbstractVector{TT}, r₀::TT, N_traj) where {TT}
    x₀ = zeros(3, N_traj)
    for xᵢ in eachcol(x₀)
        θ = 2π * rand()
        r = r₀ * rand() * [cos(θ), sin(θ), zero(eltype(X₀))]
        xᵢ .= X₀ .+ r
    end
    return x₀
end
function _initialise_fieldlines(X₀::AbstractVector{Vector{TT}}, r₀::Nothing, N_traj) where {TT}
    x₀ = hcat(X₀...)
    return x₀
end


"""
Check if the tracer has reached some ζ - ζ₀ == 0 where
  ζ = atan(u[2], u[1]) inside the event.
Taking
``arctan(y/x) - \\zeta``
Note that this is equivalent to choosing the toroidal angle.
The `mod2pi` is to shift the interval of arctan to [0,2π)
"""
function poincare_event(u, t, integrator, ζ)
    atan(u[2], u[1]) - ζ
end

function poincare_event_ζ_affect!(integrator) end



"""
Construct a `PoincarePlane` at a given ζ₀∈[0,2π) assuming there is a toroidal angle with intial points centred at `X₀` with a radius `r₀`
currently only computes a single Poincare plane
"""
function construct_poincare(coilset::CoilSet{TT}, X₀, r₀; ζ₀=zero(TT), N_traj=100, t_f=800, initial_region=nothing) where {TT}

    # We will initialise about the point X₀
    x₀ = _initialise_fieldlines(X₀, r₀, N_traj)

    if isnothing(r₀)
        N_traj = length(X₀)
    end

    # Function for moving to the next initial condition
    function prob_fn(prob, i, repeat)
        remake(prob, u0=x₀[:, i])
    end


    cb = ContinuousCallback((u, t, ∫) -> poincare_event(u, t, ∫, ζ₀), poincare_event_ζ_affect!, save_positions=(true, false))
    # Construct the problem and solve the trajectories in parallel
    ζ = (0.0, t_f / 2)
    P = ODEProblem((ẋ, x, p, t) -> field_line!(ẋ, x, p, t, coilset), x₀[:, 1], ζ)
    EP = EnsembleProblem(P, prob_func=prob_fn)
    simf = solve(EP, Tsit5(), EnsembleThreads(), trajectories=N_traj, reltol=1e-10, callback=cb, save_everystep=false, save_start=false, save_end=false)

    cb = ContinuousCallback((u, t, ∫) -> poincare_event(u, t, ∫, ζ₀), poincare_event_ζ_affect!, save_positions=(true, false))
    ζ = (0.0, -t_f / 2)
    P = ODEProblem((ẋ, x, p, t) -> field_line!(ẋ, x, p, t, coilset), x₀[:, 1], ζ)
    EP = EnsembleProblem(P, prob_func=prob_fn)
    simb = solve(EP, Tsit5(), EnsembleThreads(), trajectories=N_traj, reltol=1e-10, save_everystep=false, save_start=false, save_end=false)

    # Need to loop though outputs and store plane intersecetions
    # we do not know how many plane intersections we have a-priori
    # n_pts = mapreduce(length, +, simf.u)
    # n_pts += mapreduce(length, +, simb.u)
    # @show n_pts
    data = Vector{TT}[]
    # data = Vector{Vector{TT}}(zeros(n_pts, 3))
    # data = [zeros(3) for _ in 1:n_pts]
    for sim in simf.u
        for u in sim.u
            push!(data, u)
        end
    end
    for sim in simb.u
        for u in sim.u
            push!(data, u)
        end
    end
    poincare_data = PoincarePlane(data, ζ₀)

    return poincare_data
end








"""
Currently we define the intersection plane is hardcoded as `[0.0,1.0,0.0]`

TODO:
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
"""
Find the point ``X`` such that ``X - x=0`` where ``x`` is determined by following a field line for one turn
"""
function axis_diff(X, coilset, callback)
    fieldline = ODEProblem((ẋ, x, p, t) -> CoilFields.field_line!(ẋ, x, p, t, coilset), [X[1], 0.0, X[2]], (0.0, 800))
    sol = solve(fieldline, callback=callback, save_everystep=false)
    x = sol.u[end]
    return [X[1] - x[1], X[2] - x[3]]
end

"""
Find the axis of a coilset.
"""
function find_axis(X, coilset)
    event_affect!(integrator) = terminate!(integrator)
    affect_neg!(integrator) = nothing
    cb = ContinuousCallback(event, event_affect!; (affect_neg!)=nothing, save_positions=(true, false))

    nlp = NonlinearProblem((X, p) -> axis_diff(X, coilset, cb), X)
    nls = solve(nlp)
    return [nls.u[1], 0.0, nls.u[2]]
end
