
abstract type ExecutionMode end
struct Serial <: ExecutionMode end
struct Threaded <: ExecutionMode end

"""
Used to determine which method to use to evaluate ``\\mathbf{B}`` or ``\\mathbf{A}``.

Currently implemented is:
    - `CompactLinear` for [`biot_savart_compactlinearsegment`](@ref)
"""
abstract type EvaluationMode end
struct CompactLinear <: EvaluationMode end
struct Curviature <: EvaluationMode end

# struct CompactLinear end


"""
Evaluate the Biot Savart integral using the `CompactLinear()` segments.

`` B = \\frac{\\mu_0 I}{4\\pi} \\mathbf{R}_{i} \\times \\mathbf{R}_{i+1} \\frac{R_i R_{i+1}}{R_iR_{i+1} (R_iR_{i+1} + \\mathbf{R}_i\\cdot\\mathbf{R}_{i+1})} ``

The expression segments from the [ABCXYZ code](https://www.osti.gov/biblio/7130165) (equation 31), referenced in [Hanson and Hirshman (2002)](https://pubs.aip.org/pop/article/9/10/4410/265182/Compact-expressions-for-the-Biot-Savart-fields-of)
"""
function biot_savart_compactlinearsegment(pt1, pt2, X)
    Rᵢ = pt1 - X
    Rᵢ₊₁ = pt2 - X
    rᵢ = norm(Rᵢ)
    rᵢ₊₁ = norm(Rᵢ₊₁)
    RcrossR = cross(Rᵢ, Rᵢ₊₁)
    B = RcrossR * (rᵢ + rᵢ₊₁) / (rᵢ * rᵢ₊₁ * (rᵢ * rᵢ₊₁ + dot(Rᵢ, Rᵢ₊₁)))
    return B
end



"""
    biot_savart(c, X, evaluation_mode)

Evaluates the Biot Savart integral to compute ``\\mathbf{B}`` using a `Coil` or `AbstractCoilSet`, a position or set of positions `X` and
using some [`EvaluationMode`](@ref)
"""
function biot_savart end

# Single point single coil
function biot_savart(coil::Coil{TT,GEO}, X::Vector{TT}, ::CompactLinear) where {TT<:Real,GEO}
    B = mapreduce(i -> biot_savart_compactlinearsegment(coil[i], coil[i+1], X), +, 1:coil.length-1)
    Bfac = (coil.J * μ₀ / 4π)::TT
    Bscal = B * Bfac
    return Bscal
end

# Many coils one point
biot_savart(coilset::CoilSet, X::Vector{TT}, evaluation_mode::EvaluationMode) where {TT<:Real} =
    mapreduce(coil -> biot_savart(coil, X, evaluation_mode), +, coilset)

# Composite coil set - should lower to the one above -- output may not be type stable??
biot_savart(ccs::CompositeCoilSet, X::Vector{TT}, evaluation_mode) where {TT<:Real} =
    mapfoldl(cset -> biot_savart(cset, X, evaluation_mode), +, ccs.Group)

# Many coils many points
biot_savart(coilset::AbstractCoilSet, X::AbstractArray{Vector{TT}}, evaluation_mode::EvaluationMode) where {TT} =
    map(pt -> biot_savart(coilset, pt, evaluation_mode), X)




"""
Evaluate the Biot Savart integral for the vector potential ``A`` using the `CompactLinear` segments from Hanson and Hirshman 2002, equation 7.
"""
function biot_savart_A!(A, coil::Coil{TT,GEO}, X, ::CompactLinear) where {TT,GEO<:Tuple}
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
# function biot_savart!(B::Vector{TT}, coil::Coil{TT,PointCurvature}, X::Vector{TT}, ::Curviature) where {TT}

#     # α = κ * abs(δr′)^2 / 12

#     for I in 1:coil.length-1
#         Rᵢ = coil.Geometry.Geometry[I] - X
#         Rᵢ₊₁ = coil.Geometry.Geometry[I+1] - X
#     end
#     B .*= coil.J * μ₀ / 4π
# end






"""
In place calls
"""
function biot_savart! end

# Multiple coils single point in place
function biot_savart!(B::Vector{TT}, coilset::AbstractCoilSet, X::Vector{TT}, evaluation_mode, α=zero(TT)) where {TT<:Real}
    @. B = α * B
    B += biot_savart(coilset, X, evaluation_mode)
end
# Multiple coils multiple points
function biot_savart!(B::AbstractArray{Vector{TT}}, coilset::AbstractCoilSet, X, evaluation_mode, α=zero(TT)) where {TT<:Real}
    @. B = α * B
    for i in eachindex(B)
        B[i] += biot_savart(coilset, X[i], evaluation_mode)
    end
end

# function biot_savart!(B::AbstractArray{Vector{TT}}, ccs::CompositeCoilSet, X, evaluation_mode, α=zero(TT)) where {TT<:Real}
#     @. B = α * B
#     for i in eachindex(B)
#         B[i] .+= biot_savart(ccs, X[i], evaluation_mode)
#     end
# end






# biot_savart!(B,coil,X,::Evaluation{:CompactLinear}) = BiotSavart_CompactLinear!()



biot_savart_A!(B::Vector{Vector{TT}}, coil, X, evaluation_mode) where {TT} = map((b, pt) -> biot_savart_A!(b, coil, pt, evaluation_mode), B, X)






biot_savart_A(coil::Coil, X::Vector{TT}, evaluation_mode) where {TT} = biot_savart_A!(zeros(3), coil, X, evaluation_mode)
biot_savart_A(coil::Coil, X::Vector{Vector{TT}}, evaluation_mode) where {TT} = map!(pt -> biot_savart_A(coil, pt, evaluation_mode), [zeros(3) for _ in eachindex(X)], X)
