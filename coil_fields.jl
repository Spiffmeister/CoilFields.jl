module coil_fields



using SpecialFunctions: ellipk, ellipe
using LinearAlgebra
using Contour



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

    dkdr = k / 2r - (r-r′) * k / ((r + r′)^2 + (z - z′)^2)

    dGdr = G(r,z,r′,z′) / 2r + dkdr * sqrt((r + r′)^2 + (z - z′)^2) / (4π*k) * (-2*ellipk(1-tmpvar) + (2-k²)/(1-k²) * ellipe(k²) )

    return dGdr
end

function dGdz(r,z,r′,z′)
    tmpvar = one_mk_sqr(r,z,r′,z′)
    k² = 1 - tmpvar
    k = sqrt(k²)

    dkdz = -(z-z′)*k / ((r + r′)^2 + (z - z′)^2)

    dGdz = dkdz * sqrt((r + r′)^2 + (z - z′)^2)/(4π*k) * (-2*ellipk(1-tmpvar) + (2-k²)/(1-k²) * ellipe(tmpvar))

    return dGdz
end






ψ(coil::Coil,R,Z) = G(R,Z,coil.R,coil.Z) * μ₀ * coil.J


function ψ(coils::Array{Coil{TT}}, R, Z) where TT
    field_value = TT(0)
    # @show R, Z
    for coil in coils
        # @show field_value
        field_value += ψ(coil,R,Z)
        # @show field_value
    end
    return field_value
end




dψdr(coil,R,Z) = dGdr(R,Z,coil.R,coil.Z) * μ₀ * coil.J

dψdz(coil,R,Z) = dGdz(R,Z,coil.R,coil.Z) * μ₀ * coil.J



∇ψ(coil,R,Z) = dψdr(coil,R,Z), dψdz(coil,R,Z)

function ∇ψ(coils::Array{Coil{TT}},R,Z) where TT
    dpdr = zero(TT)
    dpdz = zero(TT)

    for coil in coils
        dpdr += dψdr(coil,R,Z)
        dpdz += dψdr(coil,R,Z)
    end

    return dpdr, dpdz
end









#= MAGNETIC FIELD FUNCTIONS =#


"""
Compute the poloidal magnetic field from a group of coils
"""
function B_poloidal(coils::Array{Coil{TT}},R,Z) where TT
    dpdr, dpdz = ∇ψ(coils,R,Z)

    return -dpdz/R, dpdr/R
end









#= FLUX SURFACE FINDING =#

# function find_flux_surface!(R,Z,field,coils;xlims=[-5.0,5.0],ylims=[-5.0,5.0])
# end
"""
Finds a flux surface with a given value of ψ
"""
function find_flux_surface(R::TT,Z::TT,gridx,gridy,field::Matrix{TT},coils::Array{Coil{TT}}) where TT

    local flux_surface

    contour_value = coil_fields.ψ(coils,R,Z)
    contour_trace = Contour.contour(gridx,gridy,field,contour_value)


    # We want the periodic contour
    for line in contour_trace.lines
        if line.vertices[1] == line.vertices[end]
            flux_surface = line
            continue
        end
    end

    flux_surface_verts = hcat([[X[1],X[2]] for X in flux_surface.vertices]...)

    segment_lengths = [norm(flux_surface_verts[:,i+1] - flux_surface_verts[:,i]) for i in Base.OneTo(length(flux_surface.vertices)-1)]

    curve_length = sum(segment_lengths)

    t(i) = segment_lengths[i]/curve_length




end


struct Fourier{TT,MODES}
    R_cosine    ::Vector{TT}
    R_sine      ::Vector{TT}
    Z_cosine    ::Vector{TT}
    Z_sine      ::Vector{TT}
    Fourier(m) = new{Float64,m}(zeros(Float64,m+1),zeros(Float64,m+1),zeros(Float64,m+1),zeros(Float64,m+1))
end

function fourier_flux_surface(flux_surface, number_of_modes, symmetric=true)

    FSeries = Fourier(number_of_modes)

    xlims = minimum(flux_surface[1,:]), maximum(flux_surface[1,:])
    ylims = minimum(flux_surface[2,:]), maximum(flux_surface[2,:])
    
    FSeries.R_cosine[1] = +(xlims...) * 0.5
    FSeries.R_cosine[2] = -(xlims...) * 0.4 #TODO: why do this scaling?
    FSeries.Z_sine[2] = ylims[2]

    for mode in Base.OneTo(number_of_modes)
        FSeries.R_cosine
        if symmetric
        end
    end

end






end