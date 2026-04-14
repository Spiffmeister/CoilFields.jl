"""
"""
struct PointCurvature{TT,VT} <: AbstractCoilGeometry
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
