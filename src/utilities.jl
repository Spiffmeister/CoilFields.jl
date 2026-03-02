point_axis(coil::Coil{TT,<:NTuple},axis) where {TT} = (pt[axis] for pt in coil.Geometry)

"""
Determine the `centre_of_mass` of the coil, given by the average in each Cartesian direction.
"""
centre_of_mass(coil::Coil{TT,<:NTuple}) where {TT} = (reduce(+, point_axis(coil, AX))/TT(coil.length) for AX in 1:3) |> collect



"""
Find a plane which lies in the "average" plane of the coil

Computed by finding the `centre_of_mass` to use as an origin and then computing ``n = v_1 \times v_2`` for each point.
"""
function average_plane(coil::Coil)

    com = centre_of_mass(coil)
    # Skip every length/10 points to make sure the triangles are not too thin
    skip = fld(coil.length, 10)

    # Compute a rolling average of normalised vectors
    n = zeros(eltype(coil), 3)
    for i in eachindex(coil)
        v₁ = com - coil[i]
        v₂ = com - coil[mod1(i + skip, coil.length)]
        nᵢ = cross(v₁, v₂)
        n .+= nᵢ / norm(nᵢ)
    end
    n ./= coil.length

    return n
end
