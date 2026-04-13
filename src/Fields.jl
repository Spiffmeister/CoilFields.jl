
abstract type ExecutionMode end
struct Serial <: ExecutionMode end
struct Threaded <: ExecutionMode end


# struct Evaluation{EVAL_TYPE} end
# const Linear = Evaluation{:LinearSegment}() ## TODO: Implement

abstract type EvaluationMode end
struct CompactLinear <: EvaluationMode end
struct Curviature <: EvaluationMode end

# struct CompactLinear end


"""
CompactLinear segments from [Hanson and hirshman (2002)](https://pubs.aip.org/pop/article/9/10/4410/265182/Compact-expressions-for-the-Biot-Savart-fields-of)
"""
# const CompactLinear = Evaluation{:CompactLinear}()
# const LinearCurviature = Evaluation{:LinearCurviature}() ## TODO: Implement

"""
Evaluate the Biot Savart integral using the `CompactLinear` segments from

Each segment is computed using the analytic form of the Biot Savart integral,
``\\int_0^1``
"""
function Biot_Savart_CompactLinearSegment(pt1, pt2, X)
    Rᵢ = pt1 - X
    Rᵢ₊₁ = pt2 - X
    rᵢ = norm(Rᵢ)
    rᵢ₊₁ = norm(Rᵢ₊₁)
    RcrossR = cross(Rᵢ, Rᵢ₊₁)
    B = RcrossR * (rᵢ + rᵢ₊₁) / (rᵢ * rᵢ₊₁ * (rᵢ * rᵢ₊₁ + dot(Rᵢ, Rᵢ₊₁)))
    return B
end
# function Biot_Savart!(B::Vector{TT}, coil::Coil, X::Vector{TT}, ::CompactLinear) where {TT<:Real}
#     for i in 1:coil.length-1
#         Biot_Savart_CompactLinearSegment!(B, coil[i], coil[i+1], X)
#     end
#     B .*= (coil.J * μ₀ / 4π)
#     return B
# end



"""
Biot_Savart
"""
function Biot_Savart end

# Single point single coil
function Biot_Savart(coil::Coil{TT,GEO}, X::Vector{TT}, ::CompactLinear) where {TT<:Real,GEO}
    B = mapreduce(i -> Biot_Savart_CompactLinearSegment(coil[i], coil[i+1], X), +, 1:coil.length-1)
    B = B * (coil.J * μ₀ / 4π)
    return Tuple(B)
end

"""
Evaluate the Biot Savart integral for the vector potential ``A`` using the `CompactLinear` segments from Hanson and Hirshman 2002, equation 7.
"""
function Biot_Savart_A!(A, coil::Coil{TT,GEO}, X, ::CompactLinear) where {TT,GEO<:Tuple}
    for I in 1:coil.length-1
        xᵢ = coil.Geometry[I]
        xᵢ₊₁ = coil.Geometry[I+1]

        Rᵢ = X - xᵢ
        Rᵢ₊₁ = X - xᵢ₊₁
        L = norm(xᵢ₊₁ - xᵢ)

        ê = (Rᵢ₊₁ - Rᵢ) / L

        ϵ = L / (norm(xᵢ) + norm(xᵢ₊₁))
        A .+= ê * log((1 + ϵ) / (1 - ϵ))
    end
    A .*= coil.J * μ₀ / 4π
end




"""
Evaluate the Biot Savart integral using the second order curvature form from
"""
function Biot_Savart!(B::Vector{TT}, coil::Coil{TT,PointCurvature}, X::Vector{TT}, ::Curviature) where {TT}

    # α = κ * abs(δr′)^2 / 12

    for I in 1:coil.length-1
        Rᵢ = coil.Geometry.Geometry[I] - X
        Rᵢ₊₁ = coil.Geometry.Geometry[I+1] - X
    end
    B .*= coil.J * μ₀ / 4π
end



# Many coils one point
Biot_Savart(coilset::CoilSet, X::Vector{TT}, evaluation_mode::EvaluationMode) where {TT<:Real} =
    mapreduce(coil -> Biot_Savart(coil, X, evaluation_mode), .+, coilset)
# Composite coil set - should lower to the one above
Biot_Savart(ccs::CompositeCoilSet, X::Vector{TT}, evaluation_mode) where {TT<:Real} =
    mapreduce(cset -> Biot_Savart(cset, X, evaluation_mode), .+, ccs.Group)

# Many coils many points
Biot_Savart(coilset::AbstractCoilSet, X::AbstractArray{Vector{TT}}, evaluation_mode::EvaluationMode) where {TT} =
    map(pt -> Biot_Savart(coilset, pt, evaluation_mode), X)



"""
In place calls
"""
function Biot_Savart! end

# Multiple coils single point in place
function Biot_Savart!(B::Vector{TT}, coilset::AbstractCoilSet, X::Vector{TT}, evaluation_mode, α=zero(TT)) where {TT<:Real}
    @. B = α * B
    B .+= Biot_Savart(coilset, X, evaluation_mode)
end
# Multiple coils multiple points
function Biot_Savart!(B::AbstractArray{Vector{TT}}, coilset::AbstractCoilSet, X, evaluation_mode, α=zero(TT)) where {TT<:Real}
    @. B = α * B
    for i in eachindex(B)
        B[i] .+= Biot_Savart(coilset, X[i], evaluation_mode)
    end
end

# function Biot_Savart!(B::AbstractArray{Vector{TT}}, ccs::CompositeCoilSet, X, evaluation_mode, α=zero(TT)) where {TT<:Real}
#     @. B = α * B
#     for i in eachindex(B)
#         B[i] .+= Biot_Savart(ccs, X[i], evaluation_mode)
#     end
# end






# Biot_Savart!(B,coil,X,::Evaluation{:CompactLinear}) = BiotSavart_CompactLinear!()



Biot_Savart_A!(B::Vector{Vector{TT}}, coil, X, evaluation_mode) where {TT} = map((b, pt) -> Biot_Savart_A!(b, coil, pt, evaluation_mode), B, X)






Biot_Savart_A(coil::Coil, X::Vector{TT}, evaluation_mode) where {TT} = Biot_Savart_A!(zeros(3), coil, X, evaluation_mode)
Biot_Savart_A(coil::Coil, X::Vector{Vector{TT}}, evaluation_mode) where {TT} = map!(pt -> Biot_Savart_A(coil, pt, evaluation_mode), [zeros(3) for _ in eachindex(X)], X)
