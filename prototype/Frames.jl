module Frames

using LinearAlgebra: norm, dot


@doc raw"""
Object for the G-frame used to construct a GVEC G-frame.

```julia
f = GFrame()
```

The curve can be evaluated at any point ``\zeta\in[0,2\pi)``
```math
f(q) = \mathbf{X}_0(q^3) + q^1 \mathbf{N}(q^3) + q^2 \sigma \mathbf{B}(q^3)
```
"""
struct GFrame{XCURVE,NCURVE,BCURVE}
    X₀::XCURVE
    N::NCURVE
    B::BCURVE

    function GFrame(X, N, B)
        return new{typeof(X),typeof(N),typeof(B)}(X, N, B)
    end
end


"""
Apply the frame mapping
``h : (q^1, q^2, q^3) \\mapsto (x,y,z)``
"""
function apply_frame(q::AbstractVector, frame::GFrame)
    return frame.X₀(q[3]) + q[1] * frame.N(q[3]) + q[2] * frame.B(q[3])
end





@doc raw"""
Returns the metric at a point ``q=(q^1, q^2, q^3)``

The metric is given by
```math
\begin{align}
G = (Dh)^T Dh =
\begin{bmatrix}
|\mathbf{N}|^2 & \mathbf{N}\cdot\mathbf{B} & \mathbf{N}\cdot\tilde{\mathbf{T}}_q \\
G_{12} & |\mathbf{B}|^2 & \mathbf{B}\cdot\tilde\mathbf{T}}_q \\
G_{13} & G_{23} & |\tilde{\mathbf{T}}_q|^2
\end{bmatrix}
\end{align}
```
"""
function metric(q::AbstractVector, frame::GFrame)
    g = zeros(eltype{q}, 3, 3) #allocating, remove

    # T = frame.T(q)
    T_q, N, B = ∂x_∂qⁱ(q[3], frame)
    # N = frame.N(q)
    # B = frame.B(q)

    g[1, 1] = norm(N)^2
    g[1, 2] = g[2, 1] = dot(N, B)
    g[1, 3] = g[3, 1] = dot(N, T_q)

    g[2, 2] = norm(B)^2
    g[2, 3] = g[3, 2] = dot(B, T_q)

    g[3, 3] = norm(T_q)^2

    return g
end

@doc raw"""
Compute the metric elements

The curves N and B must have derivatives implemented which can be evaluated at a point.

```julia
derivative(N, ζ)
```
is
```math
\left. \frac{\partial N}{\partial \zeta} \right|_{q_3}
```

"""
function ∂x_∂qⁱ(q, frame::GFrame)
    ζ = q[3]
    ∂x∂q¹ = frame.N(ζ)
    ∂x∂q² = frame.B(ζ)
    ∂x∂q³ = frame.X₀(ζ) + derivative(frame.N, ζ) + derivative(frame.B, ζ)
    return ∂x∂q¹, ∂x∂q², ∂x∂q³
end




end
