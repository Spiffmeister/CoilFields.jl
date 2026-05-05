using CoilFields
using Frames


"""
Evaluate ``\\mathbf{B} \\cdot \\mathbf{n}`` at a point `q=(q¹,q²,q³)`

The `surface` object requires a `surface_normal` function to be implemented.
"""
function Bdotn(q, coilset, surface, frame)
    n = surface_normal(q, surface, frame)
    x = frame_mapping(q, frame)
    B = biot_savart(coilset, x, CompactLinear())

    return dot(B, n)
end


"""
Evaluate the quadratic flux
```math
\int_{\partial \Omega} \frac{\|\mathbf{B}\cdot\bm{n}\|^2 - B_T}{\sigma} \text{d} S
```
where ``\\sigma`` is a normalisation factor and ``B_T`` is a target value for the flux.

TODO: Implement `σ` and `B_T`.
"""
function quadratic_flux(q::AbstractVector{AbstractVector}, coilset::CoilFields.AbstractCoilSet, surface, frame::Frames.GFrame)
    qf = zero(eltype(coilset))
    for qᵢ in q
        qf += Bdotn(qᵢ, coilset, surface, frame)^2
    end
end
