
# abstract type Coil end


const μ₀ = 4π * 1e-7


struct Coil{TT,GEOMETRY}
    Geometry    :: GEOMETRY
    J           :: TT
end

struct CoilSet
    Coils :: Vector{Coil}
end







struct evaluation{EVAL_TYPE} end
const Linear = evaluation{:LinearSegment}()




"""
Evaluate the Biot Savart integral using linear segments

Each segment is computed using the analytic form of the Biot Savart integral,
``\\int_0^1``
"""
function Biot_Savart(coil{TT,GEO},X::Vector{TT},::evaluation{:Linear}) where {TT,GEO<:AbstractVector{AbstractVector{TT}}}
    Bx = By = Bz = TT(0)
    for I in 1:coil.length-1
        
    end
    return coil.J * μ₀ * (Bx, By, Bz)
end






