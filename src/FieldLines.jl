
struct FieldLine{TT}
    position::Vector{Vector{TT}}
    t::Vector{TT}
end

struct PoincarePlane{TT}
    points::Vector{Vector{TT}}
    ζ::TT
end

function field_line!(ẋ, x, p, t, coilset)
    B = Biot_Savart(coilset, x, CompactLinear())
    ẋ .= B
    ẋ ./= norm(B)
end
function field_line_RZ!(ẋ, x, p, t, coilset)
    X = (x[1] * cos(t), x[1] * sin(t), x[2])
    B = Biot_Savart(coilset, X, CompactLinear())
    R = x[1]
    Z = x[2]
    e_ρ = (cos(t), sin(t), 0)
    e_ϕ = (-sin(t), cos(t), 0)

    B_ρ = dot(B, e_ρ)
    B_ϕ = dot(B, e_ϕ)

    ẋ .= (B_ρ, B[3])
    ẋ .*= R / B_ϕ
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
    x₀ = zeros(length(X₀), N_traj)
    for xᵢ in eachcol(x₀)
        θ = 2π * rand()
        xᵢ[1] = X₀[1] .+ r₀ * rand() * cos(θ)
        xᵢ[2] = X₀[2] .+ r₀ * rand() * sin(θ)
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




# Function for moving to the next initial condition
function _prob_fn(prob, i, repeat, x₀)
    remake(prob, u0=x₀[:, i])
end
# Actual function call when using the field line tracing
function _field_line_trace(x₀, ζ, N_traj, coilset, ζ₀, integrator, events=nothing, save_times=[])

    if !isnothing(events)
        cb = ContinuousCallback((u, t, ∫) -> poincare_event(u, t, ∫, ζ₀), poincare_event_ζ_affect!, save_positions=(true, false))
    else
        cb = nothing
    end

    P = ODEProblem((ẋ, x, p, t) -> field_line!(ẋ, x, p, t, coilset), x₀[:, 1], ζ)
    EP = EnsembleProblem(P, prob_func=(prob, i, repeat) -> _prob_fn(prob, i, repeat, x₀))
    sim = solve(EP, integrator, EnsembleThreads(), trajectories=N_traj, reltol=1e-10, callback=cb, save_everystep=false, save_start=false, save_end=false, saveat=save_times)
    return sim
end


"""
Construct a `PoincarePlane` at a given ζ₀∈[0,2π) assuming there is a toroidal angle with intial points centred at `X₀` with a radius `r₀`
currently only computes a single Poincare plane
"""
function construct_poincare(coilset, X₀, r₀; ζ₀=zero(eltype(coilset)), event=nothing, saveat=[], N_traj=100, t_f=800, integrator=Tsit5())

    # We will initialise about the point X₀
    x₀ = _initialise_fieldlines(X₀, r₀, N_traj)

    if isnothing(r₀)
        N_traj = length(X₀)
    end

    # Function for moving to the next initial condition
    function prob_fn(prob, i, repeat)
        remake(prob, u0=x₀[:, i])
    end

    if isnothing(event) && isempty(saveat)
        saveat = collect(0:2π:t_f)
    end

    # Split the trajectories into forward and backward tracing
    simf = _field_line_trace(x₀, (0, t_f / 2), N_traj, coilset, ζ₀, integrator, event, saveat)
    simb = _field_line_trace(x₀, (0, -t_f / 2), N_traj, coilset, ζ₀, integrator, event, -saveat)


    data = Vector{eltype(coilset)}[]

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

At the moment implemented so that the initial point should be on the x-z plane
"""
function find_axis(X, coilset)
    event_affect!(integrator) = terminate!(integrator)
    affect_neg!(integrator) = nothing
    cb = ContinuousCallback(event, event_affect!; (affect_neg!)=nothing, save_positions=(true, false))

    nlp = NonlinearProblem((X, p) -> axis_diff(X, coilset, cb), X)
    nls = solve(nlp)
    return [nls.u[1], 0.0, nls.u[2]]
end
