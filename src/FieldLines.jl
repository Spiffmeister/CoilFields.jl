
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


"""
Check if the tracer has reached some ζ - ζ₀ == 0 where
  ζ = atan(u[2], u[1]) inside the event.
Note that this is equivalent to choosing the toroidal angle.
The `mod2pi` is to shift the interval of arctan to [0,2π)
"""
function poincare_event(u, t, integrator, ζ)
    atan(u[2], u[1]) - ζ
end

function poincare_event_ζ_affect!(integrator) end


"""
Construct a `PoincarePlane` at a given ζ₀∈[0,2π) assuming there is a toroidal angle
"""
function construct_poincare(coilset::CoilSet{TT}; ζ₀=zero(TT), N_traj=100) where {TT}

    # cb =

end
