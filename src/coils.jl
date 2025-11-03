
# abstract type Coil end
# abstract type CoilGeometry end

const μ₀ = 4π * 1e-7

"""
    Coil{TT,GEOMETRY}

Geometry: Coil geometry
J: current
length: Coil resolution
"""
struct Coil{TT,GEOMETRY}
    Geometry::GEOMETRY
    J::TT
    length::Int
end


"""
"""
struct PointCurvature{TT,VT}
    points::VT
    normal::VT
    κ::Vector{TT}
end
# function PointCurvature(ℱ::Fourier{TT}, θ, ζ)
#     xyz = []
#     κ = []
#     new{TT,typeof(xyz)}(xyz, κ)
# end

struct ToroidalCoil{VT<:AbstractVector}
    point::VT
end

"""
    CoilSet
Store coils as a vector of object [`Coil`](@ref)

TODO: Make GEOMETRY flexible
"""
struct CoilSet{TT,GEOMETRY}
    Coils::Vector{Coil{TT,GEOMETRY}}
end
Base.iterate(CS::CoilSet, state=1) = state > length(CS.Coils) ? nothing : (CS.Coils[state], state + 1)
