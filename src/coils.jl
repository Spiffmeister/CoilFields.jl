
# abstract type Coil end


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


"""
    CoilSet
Store coils as a vector of object `Coil`
"""
struct CoilSet{TT,GEOMETRY}
    Coils::Vector{Coil{TT,GEOMETRY}}
end
Base.iterate(CS::CoilSet, state=1) = state > length(CS.Coils) ? nothing : (CS.Coils[state], state + 1)







struct evaluation{EVAL_TYPE} end
# const Linear = evaluation{:LinearSegment}() ## TODO: Implement
const CompactLinear = evaluation{:CompactLinear}()
# const LinearCurviature = evaluation{:LinearCurviature}() ## TODO: Implement




"""
Evaluate the Biot Savart integral using the `CompactLinear` segments from

Each segment is computed using the analytic form of the Biot Savart integral,
``\\int_0^1``
"""
function Biot_Savart!(B::Vector{TT}, coil::Coil{TT,GEO}, X::Vector{TT}, ::evaluation{:CompactLinear}) where {TT,GEO<:AbstractVector{<:AbstractVector{TT}}}
    for I in 1:coil.length-1
        Rᵢ = coil.Geometry[I] - X
        Rᵢ₊₁ = coil.Geometry[I+1] - X

        rᵢ = norm(Rᵢ)
        rᵢ₊₁ = norm(Rᵢ₊₁)

        B .+= cross(Rᵢ, Rᵢ₊₁) * (rᵢ + rᵢ₊₁) / (rᵢ * rᵢ₊₁ * (rᵢ * rᵢ₊₁ + dot(Rᵢ, Rᵢ₊₁)))
    end
    B .*= (coil.J * μ₀ / 4π)
end

"""
Evaluate the Biot Savart integral for the vector potential ``A`` using the `CompactLinear` segments from Hanson and Hirshman 2002, equation 7.
"""
function Biot_Savart_A!(A, coil::Coil{TT,GEO}, X, ::evaluation{:CompactLinear}) where {TT,GEO<:AbstractVector{<:AbstractVector{TT}}}
    for I in 1:coil.length-1
        xᵢ = coil.Geometry[I]
        xᵢ₊₁ = coil.Geometry[I+1]

        L = norm(xᵢ₊₁ - xᵢ)

        ê = (xᵢ - xᵢ₊₁) / L

        ϵ = L / (norm(xᵢ) + norm(xᵢ₊₁))
        A .+= ê * log((1 + ϵ) / (1 - ϵ))
    end
    A .*= coil.J * μ₀ / 4π
end




"""
Evaluate the Biot Savart integral using the second order curvature form from
"""
function Biot_Savart!(B::Vector{TT}, coil::Coil{TT,PointCurvature}, X::Vector{TT}, ::evaluation{:Curviature}) where {TT}

    # α = κ * abs(δr′)^2 / 12

    for I in 1:coil.length-1
        Rᵢ = coil.Geometry.Geometry[I] - X
        Rᵢ₊₁ = coil.Geometry.Geometry[I+1] - X
    end
    B .*= coil.J * μ₀ / 4π
end


function Biot_Savart!(B::Vector{TT}, coilset::CoilSet, X, evaluation_mode) where {TT<:Real}
    tmp = zeros(eltype(B), 3)
    for coil in coilset
        B .+= Biot_Savart!(tmp, coil, X, evaluation_mode)
    end
end
function Biot_Savart(coilset::CoilSet, X, evaluation_mode)
    tmp = zeros(eltype(X), 3)
    B = zeros(eltype(X), 3)
    for coil in coilset
        B .+= Biot_Savart!(tmp, coil, X, evaluation_mode)
    end
    return B
end


# Biot_Savart!(B,coil,X,::evaluation{:CompactLinear}) = BiotSavart_CompactLinear!()

Biot_Savart!(B::Vector{Vector{TT}}, coil, X, evaluation_mode) where {TT} = tmap((b, pt) -> Biot_Savart!(b, coil, pt, evaluation_mode), B, X)
Biot_Savart_A!(B::Vector{Vector{TT}}, coil, X, evaluation_mode) where {TT} = tmap((b, pt) -> Biot_Savart_A!(b, coil, pt, evaluation_mode), B, X)

"""
    Biot_Savart(coil, X, evaluation_mode)

Evaluate the Biot Savart integral at location(s) `X` where `X` is a Cartesian point or a vector of Cartesian points (`Vector{Vector{TT}}`)
"""
Biot_Savart(coil::Coil, X::Vector{TT}, evaluation_mode) where {TT<:Real} = Biot_Savart!(zeros(3), coil, X, evaluation_mode)
Biot_Savart(coil::Coil, X::Vector{Vector{TT}}, evaluation_mode) where {TT} = map!(pt -> Biot_Savart(coil, pt, evaluation_mode), [zeros(3) for _ in eachindex(X)], X)

Biot_Savart_A(coil::Coil, X::Vector{TT}, evaluation_mode) where {TT} = Biot_Savart_A!(zeros(3), coil, X, evaluation_mode)
Biot_Savart_A(coil::Coil, X::Vector{Vector{TT}}, evaluation_mode) where {TT} = map!(pt -> Biot_Savart_A(coil, pt, evaluation_mode), [zeros(3) for _ in eachindex(X)], X)
