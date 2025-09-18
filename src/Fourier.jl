
"""
(x,y,z) = A₀₀ + ∑ Aₘₙ cos(mθ - nϕ) + ∑ Bₘₙ sin(mθ - nϕ)
"""
struct Fourier1D{TT,MODES}
    cosine    ::Vector{TT}
    sine      ::Vector{TT}
    cosine_modes    :: Vector{Int}
    sine_modes      :: Vector{Int}
    Fourier1D(m) = new{Float64,m}(zeros(Float64,m),zeros(Float64,m-1),collect(0:m-1),collect(1:m-1))
end

struct FourierCurve{TT,MODES}
    R   :: Fourier1D{TT,MODES}
    Z   :: Fourier1D{TT,MODES}
    FourierCurve(m) = new{Float64,m}(Fourier1D(m), Fourier1D(m))
end


Fourier_cosine(ℱ,t) = mapreduce(ai -> ai[2]*cos((ai[1]-1)*t*2π), +, enumerate(ℱ.cosine))
Fourier_sine(ℱ,t) = mapreduce(ai -> ai[2]*sin(ai[1]*t*2π), +, enumerate(ℱ.sine)) #TODO: need to ensure that sine has less terms (skip the zero term)

(ℱ::Fourier1D)(t) = Fourier_cosine(ℱ,t) + Fourier_sine(ℱ,t)
(ℱc::FourierCurve)(t) = ℱc.R(t), ℱc.Z(t)




_order_odd = struct oddorder() end

"""
    derivative(ℱ::Fourier1D{TT,MODES},t,order=1)

Compute the `derivative` of a `Fourier` series object with `order`
"""
function derivative(ℱ::Fourier1D{TT,MODES},t,order=1) where {TT,MODES}
    dxdt = zero(TT)
    if isodd(order)
        cosine_derivative = sin
        sine_derivative = cos
    else
        cosine_derivative = cos
        sine_derivative = sin
    end
    for I in 2:MODES # First mode is A_0
        dxdt += (-(I-1)*2π)^order * ℱ.cosine[I] * cosine_derivative((I-1) * t * 2π)
    end
    for I in 1:MODES-1 # Sine starts from n=1, no n+1 mode
        dxdt += (-I * 2π) * ℱ.sine[I] * sine_derivative(I * t * 2π)
    end
    return dxdt
end


function cosine_derivative(Aᵢ,modes,t::TT,order) where TT
    dxdt = zero(TT)
    for m in modes
        dxdt += (-1)*(m*2π)^order * Aᵢ[I] * cos(t * m * 2π)
    end
    return dxdt
end
function sine_derivative(Aᵢ,modes,t::TT,order) where TT
    dxdt = zero(TT)
    for m in modes
        dxdt += (-1)*(m*2π)^(order+1) * Aᵢ[I] * sin(t * m * 2π)
    end
    return dxdt
end

function derivative(ℱ::Fourier1D{TT,MODES},t,order=1) where {TT,MODES}
    if iseven(order)
        dxdt = cosine_derivative(@view ℱ.cosine[2:end], @view ℱ.cosine_modes[2:end], t, order) + sine_derivative(ℱ.sine, ℱ.sine_modes, t, order)
    else
        dxdt = sine_derivative(@view ℱ.cosine[2:end], @view ℱ.cosine_modes[2:end], t, order) + cosine_derivative(ℱ.sine, ℱ.sine_modes, t, order)
    end
    return dxdt
end


derivative(ℱc::FourierCurve{TT},t,order=1) where TT = derivative(ℱc.R,t,order), derivative(ℱc.Z,t,order)



curvature(ℱ,t) = norm(derivative(ℱ, t, order=2))
normal(ℱ,t) = derivative(ℱ,t,order=1)./norm(derivative(ℱ,t,order=1))

