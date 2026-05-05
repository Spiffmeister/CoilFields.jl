abstract type AbstractFourierSeries end


"""
Fourier cosine or sine series object

`SType` is either `:cos` or `:sin` and will define the type of series to call.
"""
struct Fourier{TT,SType} <: AbstractFourierSeries
    amplitudes::Vector{TT}
    modes::Int
    Fourier(stype, amplitudes, m) = new{eltype(amplitudes),stype}(amplitudes, m)
end

Base.eltype(::Fourier{TT}) where {TT} = TT
seriestype(::Fourier{TT,ST}) where {TT,ST} = ST

const FourierCosineSeries{TT} = Fourier{TT,:cos}
const FourierSineSeries{TT} = Fourier{TT,:sin}

(ℱ::FourierCosineSeries{TT})(θ) where {TT} = mapreduce(ia -> ia[2] * cos((ia[1] - 1) * θ), +, enumerate(ℱ.amplitudes))
(ℱ::FourierSineSeries{TT})(θ) where {TT} = mapreduce(ia -> ia[2] * sin(ia[1] * θ), +, enumerate(ℱ.amplitudes))


function derivative(ℱ::FourierCosineSeries{TT}, θ) where {TT}
    return -mapreduce(ia -> ia[2] * (ia[1] - 1) * sin((ia[1] - 1) * θ), +, enumerate(ℱ.amplitudes))
end
function derivative(ℱ::FourierSineSeries{TT}, θ) where {TT}
    return mapreduce(ia -> ia[2] * ia[1] * cos(ia[1] * θ), +, enumerate(ℱ.amplitudes))
end


"""
Fourier cosine or sine series object for 2D objects.

Can be called with `ℱ(θ,ζ)`.

`ℱ.amplitudes` is stored as a matrix where `m` are the rows and `n` are the columns.

TODO: Implement m and n in a way which does not require storage maybe (can probably use a generator).
"""
struct Fourier2D{TT,STYPE} <: AbstractFourierSeries
    amplitudes::Matrix{TT} #m×n matrix
    M::Matrix{Int}
    N::Matrix{Int}
    N_fp::Int
    function Fourier2D(stype, A, m_max, n_max, N_fp=1)
        M = Matrix(transpose(repeat(0:m_max, 1, 2n_max + 1)))
        N = repeat(-n_max:n_max, 1, m_max)
        return new{eltype(A),stype}(A, M, N, N_fp)
    end
end
function Fourier2D(SType, Amn::Tuple, N_fp=1; m_max=0, n_max=0)

    # If these are zero determine based on input
    if iszero(m_max)
        m_max = maximum(getindex.(Amn,2))
    end
    if iszero(n_max)
        n_max = maximum(last.(Amn))
    end

    A = zeros(eltype(Amn[1][1]), m_max+1, n_max+1)
    for (a,m,n) in Amn
        A[m+1, n+1] = a
    end

    Fourier2D(SType, A, m_max, n_max, N_fp)
end



"""
Amplitudes are stored like `ℱ.amplitudes[i,j]` with rows≡`m` and cols≡`n`
"""
Base.getindex(ℱ::Fourier2D, I...) = ℱ.amplitudes[I]
function setindex!(ℱ::Fourier2D, val, ind)
    ℱ.amplitudes[ind] = val
end

(ℱ::Fourier2D{TT,:cos})(θ, ζ) where {TT} = mapreduce(i -> ℱ.amplitudes[i] * cos(ℱ.M[i] * θ - ℱ.N[i] * ℱ.N_fp * ζ), +, eachindex(ℱ.amplitudes))
(ℱ::Fourier2D{TT,:sin})(θ, ζ) where {TT} = mapreduce(i -> ℱ.amplitudes[i] * sin(TT(ℱ.M[i]) * θ - TT(ℱ.N[i]) * ℱ.N_fp * ζ), +, eachindex(ℱ.amplitudes))





"""
`FourierSeries` object which has both `sin` and `cos` series objects (or either can be `nothing`)
"""
struct FourierSeries{TT,ST,CT}
    sin::Fourier{TT,ST}
    cos::Fourier{TT,CT}

    function FourierSeries(series...)
        TT = eltype(series[1])

        elcheck = findfirst(s -> typeof(s) <: FourierCosineSeries, series)
        cosine = typeof(elcheck) <: Nothing ? nothing : series[elcheck]

        elcheck = findfirst(s -> typeof(s) <: FourierSineSeries, series)
        sine = typeof(elcheck) <: Nothing ? nothing : series[elcheck]
        new{TT,seriestype(sine),seriestype(cosine)}(sine, cosine)
    end
end


(ℱ::FourierSeries{TT,:sin,Nothing})(θ) where {TT} = ℱ.sin(θ)
(ℱ::FourierSeries{TT,Nothing,:cos})(θ) where {TT} = ℱ.cos(θ)
(ℱ::FourierSeries{TT,:sin,:cos})(θ) where {TT} = ℱ.sin(θ) + ℱ.cos(θ)






struct FourierCurve{TT}
    x::FourierSeries
    y::FourierSeries
    z::FourierSeries
    FourierCurve(Fx, Fy, Fz) = new{Float64}(Fx, Fy, Fz)
end

(ℱc::FourierCurve)(t) = ℱc.R(t), ℱc.Z(t)


# _order_odd = struct oddorder() end

# """
#     derivative(ℱ::Fourier1D{TT,MODES},t,order=1)

# Compute the `derivative` of a `Fourier` series object with `order`
# """
# function derivative(ℱ::Fourier1D{TT,MODES}, t, order=1) where {TT,MODES}
#     dxdt = zero(TT)
#     if isodd(order)
#         cosine_derivative = sin
#         sine_derivative = cos
#     else
#         cosine_derivative = cos
#         sine_derivative = sin
#     end
#     for I in 2:MODES # First mode is A_0
#         dxdt += (-(I - 1) * 2π)^order * ℱ.cosine[I] * cosine_derivative((I - 1) * t * 2π)
#     end
#     for I in 1:MODES-1 # Sine starts from n=1, no n+1 mode
#         dxdt += (-I * 2π) * ℱ.sine[I] * sine_derivative(I * t * 2π)
#     end
#     return dxdt
# end


# function cosine_derivative(Aᵢ, modes, t::TT, order) where {TT}
#     dxdt = zero(TT)
#     for m in modes
#         dxdt += (-1) * (m * 2π)^order * Aᵢ[I] * cos(t * m * 2π)
#     end
#     return dxdt
# end
# function sine_derivative(Aᵢ, modes, t::TT, order) where {TT}
#     dxdt = zero(TT)
#     for m in modes
#         dxdt += (-1) * (m * 2π)^(order + 1) * Aᵢ[I] * sin(t * m * 2π)
#     end
#     return dxdt
# end

# function derivative(ℱ::Fourier1D{TT,MODES}, t, order=1) where {TT,MODES}
#     if iseven(order)
#         dxdt = cosine_derivative(@view ℱ.cosine[2:end], @view ℱ.cosine_modes[2:end], t, order) + sine_derivative(ℱ.sine, ℱ.sine_modes, t, order)
#     else
#         dxdt = sine_derivative(@view ℱ.cosine[2:end], @view ℱ.cosine_modes[2:end], t, order) + cosine_derivative(ℱ.sine, ℱ.sine_modes, t, order)
#     end
#     return dxdt
# end


# derivative(ℱc::FourierCurve{TT}, t, order=1) where {TT} = derivative(ℱc.R, t, order), derivative(ℱc.Z, t, order)



# curvature(ℱ, t) = norm(derivative(ℱ, t, order=2))
# normal(ℱ, t) = derivative(ℱ, t, order=1) ./ norm(derivative(ℱ, t, order=1))
