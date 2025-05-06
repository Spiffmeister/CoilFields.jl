module coil_fields



using SpecialFunctions: ellipk, ellipe
using LinearAlgebra
using Contour
using Optim: optimize


export Coil, Field, Fourier1D, FourierCurve
export psi, grad_psi, B_poloidal
export find_flux_surface
export derivative, get_RZ, fourier_flux_surface



const μ₀ = 4π * 10^-7

struct CurrentType{CT} end
const CoilCurrent = CurrentType{:External}()
const PlasmaCurrent = CurrentType{:Internal}()


"""
Coil going into the page at position ``(r,z)`` and with current ``J``.
"""
struct Coil{TT<:Float64}
    R   :: TT
    Z   :: TT
    J   :: TT
    Coil(R::TT,Z,J) where TT = new{TT}(R,Z,J)
end

"""
Array of coil objects, function ``\\psi`` and the gradient of that function which returns ``(d\\psi/dr, d\\psi/dz)``
"""
struct Field{TT<:Float64,PT<:Function,GPT<:Function}
    Coils       :: Vector{Coil{TT}}
    ψ           :: PT
    ∇ψ          :: GPT
    Field(coilarray::Vector{Coil{TT}},psi,grad_psi) where TT = new{TT,typeof(psi),typeof(grad_psi)}(coilarray,psi,grad_psi)
end




one_mk_sqr(r,z,r′,z′) = ((r - r′)^2 + (z - z′)^2)/((r + r′)^2 + (z - z′)^2)


"""
    GreensFunctionIntegralFirstKind
"""
function G(r,z,r′,z′)
    #TODO: terrible variable names

    tmpvar = one_mk_sqr(r,z,r′,z′)
    k² = 1 - tmpvar

    G = sqrt( (r + r′)^2 + (z - z′)^2 )/4π * ((1+tmpvar) * ellipk(1-tmpvar) - 2 * ellipe(k²) )
    return G
end

"""
    GreensFunctionIntegralFirstKind_dr
Derivative of G with respect to r
"""
function dGdr(r,z,r′,z′)
    tmpvar = one_mk_sqr(r,z,r′,z′)
    k² = 1 - tmpvar

    k = sqrt(k²)

    dkdr = k / 2r - (r+r′) * k / ((r + r′)^2 + (z - z′)^2)

    dGdr = G(r,z,r′,z′) / 2r + dkdr * sqrt((r + r′)^2 + (z - z′)^2) / (4π*k) * (-2*ellipk(1-tmpvar) + (2-k²)/(1-k²) * ellipe(k²) )

    return dGdr
end

function dGdz(r,z,r′,z′)
    tmpvar = one_mk_sqr(r,z,r′,z′)
    k² = 1 - tmpvar
    k = sqrt(k²)

    dkdz = -(z-z′)*k / ((r + r′)^2 + (z - z′)^2)

    dGdz = dkdz * sqrt((r + r′)^2 + (z - z′)^2)/(4π*k) * (-2*ellipk(1-tmpvar) + (2-k²)/(1-k²) * ellipe(k²))

    return dGdz
end





"""
Compute the ``\\psi`` values given a coil or coil array
"""
psi(coil::Coil,R,Z) = G(R,Z,coil.R,coil.Z) * μ₀ * coil.J
function psi(coils::Array{Coil{TT}}, R, Z) where TT
    field_value = TT(0)
    for coil in coils
        field_value += psi(coil,R,Z)
    end
    return field_value
end
function psi(field::Field{TT,PT,GPT},R,Z) where {TT,PT,GPT}
    field_value = psi(field.Coils,R,Z)
    field_value += field.ψ(R,Z)
    return field_value
end



"""
Derivatives of the ``\\psi`` field
"""
dψdr(coil,R,Z) = dGdr(R,Z,coil.R,coil.Z) * μ₀ * coil.J

dψdz(coil,R,Z) = dGdz(R,Z,coil.R,coil.Z) * μ₀ * coil.J

grad_psi(coil,R,Z) = dψdr(coil,R,Z), dψdz(coil,R,Z)

function grad_psi(coils::Array{Coil{TT}},R,Z) where TT
    dpdr = zero(TT)
    dpdz = zero(TT)

    for coil in coils
        dpdr += dψdr(coil,R,Z)
        dpdz += dψdz(coil,R,Z)
    end

    return dpdr, dpdz
end
function grad_psi(field::Field,R,Z)
    dpdr, dpdz = grad_psi(field.Coils,R,Z)
    dpdr_tmp, dpdz_tmp = field.∇ψ(R,Z)
    dpdr += dpdr_tmp
    dpdz += dpdz_tmp
    return dpdz, dpdz
end








#= MAGNETIC FIELD FUNCTIONS =#


"""
Compute the poloidal magnetic field from a group of coils

Optional parameter can be the gradient of some ``\\psi`` function which returns the derivatives (dψ/dr,dψ/dz) at a point as a vector or tuple.
"""
function B_poloidal(coils::Array{Coil{TT}},R,Z) where {TT}
    dpdr, dpdz = grad_psi(coils,R,Z)
    return -dpdz/R, dpdr/R
end
function B_poloidal(field::Field{TT,PT,GPT},R,Z) where {TT,PT,GPT}
    B_pr, B_pz = B_poloidal(field.Coils,R,Z)
    dpdr, dpdz = field.∇ψ(R,Z)
    B_pr += -dpdz/R
    B_pz += dpdr/R
    return B_pr, B_pz
end

"""
Compute the toroidal current from a ``\\psi`` field
"""
function J_toroidal(field::Field{TT,PT,GPT},R,Z) where {TT,PT,GPT}

end









#= FLUX SURFACE FINDING =#

# function find_flux_surface!(R,Z,field,coils;xlims=[-5.0,5.0],ylims=[-5.0,5.0])
# end
"""
Finds a flux surface with the ``\\psi`` value at the location ``R,Z``.
"""
function find_flux_surface(R::TT,Z::TT,gridx,gridy,field::Matrix{TT},coils_or_field::FT) where {TT,FT<:Union{Array{Coil{TT}},Field}}

    local flux_surface
    old_dist = Inf
    dist = 0.0

    contour_value = psi(coils_or_field,R,Z)
    contour_trace = Contour.contour(gridx,gridy,field,contour_value)

    # minimum_distance = minimum(abs.(line.vertices .- ))

    
    for line in contour_trace.lines
        # We want the periodic contour
        if line.vertices[1] == line.vertices[end]
            # Make sure we pick the closest contour to the (R,Z) point given
            dist = mapreduce(x->norm(x .- (R,Z)), min, line.vertices)
            if dist < old_dist
                old_dist = dist
                flux_surface = line
            end
        end
    end

    # Store verts in column vector
    flux_surface_verts = hcat([[X[1],X[2]] for X in flux_surface.vertices]...)

    return flux_surface_verts
end



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

_get_all_coefficients(ℱ::Fourier1D) = vcat(ℱ.cosine, ℱ.sine)
_get_all_coefficients(ℱc::FourierCurve) = vcat(_get_all_coefficients(ℱc.R), _get_all_coefficients(ℱc.Z))

_get_all_modes(ℱ::Fourier1D) = vcat(ℱ.cosine_modes, ℱ.sine_modes)
_get_all_modes(ℱc::FourierCurve) = vcat(_get_all_modes(ℱc.R), _get_all_modes(ℱc.Z))

_get_nonzero_modes(ℱ::Fourier1D) = vcat(ℱ.cosine .!= 0, ℱ.sine .!= 0)

function _update_coefficients(ℱ::Fourier1D{TT,MODES},coefficients::Array{TT}) where {TT,MODES}
    ℱ.cosine .= coefficients[1:MODES]
    ℱ.sine .= coefficients[MODES+1:end]
end
function _update_coefficients(ℱc::FourierCurve{TT,MODES},coefficients::Array{TT}) where {TT,MODES}
    modes = 2MODES-1
    _update_coefficients(ℱc.R,coefficients[1:modes])
    _update_coefficients(ℱc.Z,coefficients[modes+1:end])
end
# _get_all_modes(ℱ::Fourier1D)


Fourier_cosine(ℱ,t) = mapreduce(ai -> ai[2]*cos((ai[1]-1)*t*2π), +, enumerate(ℱ.cosine))
Fourier_sine(ℱ,t) = mapreduce(ai -> ai[2]*sin(ai[1]*t*2π), +, enumerate(ℱ.sine)) #TODO: need to ensure that sine has less terms (skip the zero term)

(ℱ::Fourier1D)(t) = Fourier_cosine(ℱ,t) + Fourier_sine(ℱ,t)
(ℱc::FourierCurve)(t) = ℱc.R(t), ℱc.Z(t)

function derivative(ℱ::Fourier1D{TT,MODES},t) where {TT,MODES}
    dxdt = zero(TT)
    for I in 2:MODES
        dxdt += -(I-1) * 2π * ℱ.cosine[I] * sin((I-1) * t * 2π)
    end
    for I in 1:MODES-1
        dxdt += I * 2π * ℱ.sine[I] * cos(I * t * 2π)
    end
    return dxdt
end
derivative(ℱc::FourierCurve{TT},t) where TT = derivative(ℱc.R,t), derivative(ℱc.Z,t)



"""
Get the (R,Z) coordinates of a curve defined by Fourier series
"""
function get_RZ(ℱc::FourierCurve{TT,MODES}, N; offset=zero(TT), endpoint=false) where {TT,MODES}
    endpoint ? N_offset = 1 : N_offset = 0
    θ = range(zero(TT) + offset, one(TT)+offset, N+1)[1:N+N_offset]
    R = map(ℱc.R, θ)
    Z = map(ℱc.Z, θ)
    return R,Z
end


"""
Convert a flux surface defined by a set of points to a Fourier series representation
"""
function fourier_flux_surface(flux_surface, number_of_modes, coils_or_field::CFT) where {TT,CFT<:Union{Vector{Coil{TT}}, Field}}

    curve = FourierCurve(number_of_modes)

    ψ₀ = psi(coils_or_field, flux_surface[1,1], flux_surface[2,1])

    # Find a bounding box to initialise the flux surface
    xlims = minimum(flux_surface[1,:]), maximum(flux_surface[1,:])
    ylims = minimum(flux_surface[2,:]), maximum(flux_surface[2,:])
    
    # Initial guess for the flux surface shape
    curve.R.cosine[1] = +(xlims...) * 0.5
    curve.R.cosine[2] = (xlims[2] - xlims[1]) * 0.4 #TODO: why do this scaling?
    curve.Z.sine[1] = ylims[2] * 0.8

    for mode in 2:number_of_modes

        N_I = 1+4mode #Number of integration points

        coefficients = _get_all_coefficients(curve)

        new_coefficients = optimize((coefficients) -> _minimum_psi(coefficients, curve, coils_or_field, ψ₀, N_I), coefficients)

        _update_coefficients(curve, new_coefficients.minimizer)
    end

    return curve

end


"""
Squared normalised distance to the target ``\\psi`` value
"""
function distance_to_target(ℱc::FourierCurve, coils_or_field::CFT, target, N_I) where {TT,CFT<:Union{Vector{Coil{TT}}, Field}}
    R,Z = get_RZ(ℱc, N_I)
    ψ_curve = [psi(coils_or_field, r, z) for (r,z) in zip(R,Z)]
    return (ψ_curve .- target).^2 / target.^2
end

function Bdotn(ℱc::FourierCurve, coils_or_field::CFT, N_I) where {TT,CFT<:Union{Vector{Coil{TT}}, Field}}
    θ = range(0.0,1.0,N_I+1)[1:N_I]
    R,Z = get_RZ(ℱc, N_I)
    B_p = [B_poloidal(coils_or_field,r,z) for (r,z) in zip(R,Z)]
    B_pr = first.(B_p)
    B_pz = last.(B_p)
    tmp_derivative = [derivative(ℱc, ϑ) for ϑ in θ]
    n_z = -first.(tmp_derivative) # n_z = -dr/dt
    n_r = last.(tmp_derivative) # n_r = dz/dt

    out = @. (n_r*B_pr + n_z*B_pz)^2 / ((B_pr^2 + B_pz^2)*(n_z^2 + n_r^2))
    
    return out
end

"""
Function to minimise
"""
function _minimum_psi(coefficients,ℱc,coils_or_field,ψ_target,N_I; p=1, q=p)

    _update_coefficients(ℱc, coefficients)

    mscale = 1e-3 * dot(_get_all_modes(ℱc).^(p+q), coefficients.^2) / dot(_get_all_modes(ℱc).^p, coefficients.^2)

    Δψ = sum(distance_to_target(ℱc, coils_or_field, ψ_target, N_I))/2N_I + sum(Bdotn(ℱc, coils_or_field, N_I))/2N_I + mscale
    return Δψ
end







end