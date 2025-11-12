
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



"""
Check if the tracer has reached some ζ - ζ₀ == 0 where
  ζ = atan(u[2], u[1]) inside the event.
Taking
``|arctan(y/x) - \\zeta|``
Note that this is equivalent to choosing the toroidal angle.
The `mod2pi` is to shift the interval of arctan to [0,2π)
"""
function poincare_event(u, t, integrator, ζ)
    abs(atan(u[2], u[1]) - ζ)
end

function poincare_event_ζ_affect!(integrator) end



"""
Construct a `PoincarePlane` at a given ζ₀∈[0,2π) assuming there is a toroidal angle
currently only computes a single Poincare plane
"""
function construct_poincare(coilset::CoilSet{TT}, X₀, r₀; ζ₀=zero(TT), N_traj=100, t_f=800, initial_region=nothing) where {TT}

    # ζ_slices = 2π * collect(0:N_orbs/2) .+ ζ₀

    # x₀ = rand(2, N_traj)

    # We will initialise about the point X₀
    x₀ = zeros(3, N_traj)

    for xᵢ in eachcol(x₀)
        θ = 2π * rand()
        r = r₀ * rand() * [cos(θ), sin(θ), zero(TT)]
        xᵢ .= X₀ .+ r
    end

    # Ensure things are within domains
    # x₀[1, :] = X₀[1] + rand(-r₀[1] / 2, r₀[1] / 2) #ψ
    # x₀[2, :] = X₀[2] + rand(-r₀[2] / 2, r₀[2] / 2) #θ

    # Function for moving to the next initial condition
    function prob_fn(prob, i, repeat)
        remake(prob, u0=x₀[:, i])
    end

    cb = ContinuousCallback(poincare_event, poincare_event_ζ_affect!, save_positions=(true, false))


    # Construct the problem and solve the trajectories in parallel
    ζ = (0.0, t_f / 2)
    P = ODEProblem((ẋ, x, p, t) -> field_line!(ẋ, t, x, p, coilset), x₀[:, 1], ζ)
    EP = EnsembleProblem(P, prob_func=prob_fn)
    simf = solve(EP, Tsit5(), EnsembleThreads(), trajectories=N_traj, reltol=1e-10)

    ζ = (0.0, -t_f / 2)
    P = ODEProblem((ẋ, x, p, t) -> field_line!(ẋ, t, x, p, coilset), x₀[:, 1], ζ)
    EP = EnsembleProblem(P, prob_func=prob_fn)
    simb = solve(EP, Tsit5(), EnsembleThreads(), trajectories=N_traj, reltol=1e-10)

    # Need to loop though outputs and store plane intersecetions
    # we do not know how many plane intersections we have a-priori
    n_pts = 0
    data = Vector{TT}[]
    for sim in simf.u
        # n_pts += len(sim.t)
        for u in sim.u
            push!(data, u)
        end
    end
    for sim in simb.u
        # n_pts += len(sim.t)-1 # remove initial point
        push!(data, u[2:end])
    end
    poincare_data = PoincarePlane(data, ζ₀)

    return poincare_data
end
