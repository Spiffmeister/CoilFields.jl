

include("./coil_fields.jl")
# using Revise
# using .coil_fields

using GLMakie

using Contour

using Optim





CoilRs = [0.8678435,0.8678435,4.49624,4.49624,8.638285,8.638285,1.19984,1.19984,2.62973,2.62973,9.654945,9.654945]
CoilZs = [1.860055 ,-1.860055 , 3.69586  ,-3.69586  , 2.468375 ,-2.468375 , 3.34104  ,-3.34104  , 5.249245 ,-5.249245 , 3.295555 ,-3.295555]
CoilJs = [-293555.433002323,-293555.433002323,-17349.2247843043,-17349.2247843043, 60023.3897873202, 60023.3897873202,-356061.424859307,-356061.424859307,-65883.6778882385,-65883.6778882385,-11539.671427341,-11539.671427341]
"""
R, Z , factor,  number, name, Current
  0.0  -100.0    1
  0.0      100.0       0   1  bg_BtBz_01    3.884526409876309E+06
  0.8678435   1.860055  1  2  PFsym_a_01  -293555.433002323
  0.8678435  -1.860055  1  3  PFsym_a_01  -293555.433002323
  4.49624     3.69586   1  4  PFsym_a_02  -17349.2247843043
  4.49624    -3.69586   1  5  PFsym_a_02  -17349.2247843043
  8.638285    2.468375  1  6  PFsym_a_03   60023.3897873202
  8.638285   -2.468375  1  7  PFsym_a_03   60023.3897873202
  1.19984     3.34104   1  8  PFsym_b_01  -356061.424859307
  1.19984    -3.34104   1  9  PFsym_b_01  -356061.424859307
  2.62973     5.249245  1 10  PFsym_b_02  -65883.6778882385
  2.62973    -5.249245  1 11  PFsym_b_02  -65883.6778882385
  9.654945    3.295555  1 12  PFsym_b_03  -11539.671427341
  9.654945   -3.295555  1 13  PFsym_b_03  -11539.671427341
"""


CoilArray = [coil_fields.Coil(R,Z,J) for (R,Z,J) in zip(CoilRs,CoilZs,CoilJs)]


currentarray = CoilArray




# fsurface = coil_fields.find_flux_surface(3.15, 0.0, X, Y, Z_full, currentarray)

# number_of_modes = 3

# fourier_surface = coil_fields.fourier_flux_surface(fsurface,number_of_modes,currentarray)



# ψ_field(r,z) = 5exp(-2(25(r-3)^2 + 5(z-0.0)^2))







function ψ_field(r,z)
    # Radial ψ field with

    peak = 1.0
    
    α_r, α_z    = 100.0,50.0 # Slope strength
    R₀, Z₀      = 3.0, 0.0 # Origin
    w_r         = 5.0 # Width

    R = -sqrt( (r-R₀)^2*α_r + (z-Z₀)^2*α_z) + w_r

    return peak * (tanh(R) + 1)/2
end


function gradient_psi_field(r,z)
    peak = 1.0
    α_r, α_z = 100.0, 50.0
    R₀, Z₀ = 3.0, 0.0
    w_r = 5.0

    R = -sqrt( (r-R₀)^2*α_r + (z-Z₀)^2*α_z) + w_r
    
    dRdr = -α_r*r / R
    dRdz = -α_z*z / R

    dψdr = peak * dRdr * (1 - tanh(R)^2)
    dψdz = peak * dRdz * (1 - tanh(R)^2)

    return dψdr, dψdz

end


function bisection(fn, window, tol=1e-2)
    err = 1.0
    pt = +(window...)/2
    while err > tol
        q = fn(pt)
        if q*fn(window[1]) < 0.0
            window[1] = pt
        else
            window[2] = pt
        end
        pt = +(window...)/2
        err = abs(q - fn(pt))
    end
    return pt
end


field = coil_fields.Field(currentarray, ψ_field, gradient_psi_field)

@show coil_fields.psi(field, 3.5, 0.0)
@show coil_fields.B_poloidal(field, 3.15, 0.0)

@show bisection(R->coil_fields.psi(field,R,0.0),[3.0,4.0])






X = range(0.1,11,1001)
Y = range(-7,7,1001)

Z_coils = [coil_fields.psi(currentarray,x,y) for x in X, y in Y]

Z_full = Z_coils .+ [ψ_field(x,y) for x in X, y in Y]


# contour_xy = [3.5,0.0]
# contour_value = coil_fields.psi(field,contour_xy...)
# contour_trace = Contour.contour(X,Y,Z_full,contour_value)


fsurface = coil_fields.find_flux_surface(3.5, 0.0, X, Y, Z_full,field)
fourier_surface = coil_fields.fourier_flux_surface(fsurface, 2, field)



f = Figure(); 
axf = Axis(f[1,1])
contour!(axf,X,Y,Z_full, levels=500)
scatter!(axf,CoilRs,CoilZs,color=:red)

vlines!(axf, [3.0], color=:red, linestyle=:dash)
hlines!(axf, [0.0], color=:red, linestyle=:dash)

# contour!(axf,X,Y,Z_full, levels=[1e-8], color=:black)
# xlims!(axf2,1,4); ylims!(axf2,-3,3)

lines!(axf, fsurface[1,:], fsurface[2,:], color=:black)

axf2 = Axis(f[1,2])
lines!(axf2, ψ_field.(3.0,Y), Y)
lines!(axf2, Z_full[argmin(abs.(X.-3.0)),:], Y)
hlines!(axf2, [0.0], color=:red)

axf3 = Axis(f[2,1])
axf3_cl1 = lines!(axf3, X, ψ_field.(X,0.0), label=L"$\psi$ function")
axf3_cl2 = lines!(axf3, X, Z_full[:,501], label=L"Full $\psi$")
vlines!(axf3, [3.0], color=:red)

axislegend(axf3)

Colorbar(f[1,3],limits=(minimum(Z_full),maximum(Z_full)))



# scatter!(axf,PlasmaRs,PlasmaZs,color=:black)

# scatter!(axf,contour_verts[1,:], contour_verts[2,:],color=:red)

# contour!(axf,X,Y,Z_full,levels=[coil_fields.psi(currentarray, 3.15,0.0)],color=:black)

fcurveR,fcurveZ = coil_fields.get_RZ(fourier_surface,40, endpoint=true)

scatterlines!(axf, fcurveR, fcurveZ, color=:red)


# xlims!(axf,2,5)
# ylims!(axf,-2,2)


colsize!(f.layout, 1, Relative(0.8))
rowsize!(f.layout, 1, Relative(0.8))

# axf2 = Axis(f[2,1])
# lines!(axf2, X, ψ_field.(X,0.0), label=L"\psi")
# axislegend(axf2)

# axf3 = Axis(f[1,2])
# lines!(axf3, Y, ψ_field.(3.0,Y), label=L"\psi")



f












#=

plascoil_position = Observable([Point2f(PlasmaRs[1],PlasmaZs[1])])
# plascoil2 = Observable(Point2f(PlasmaRs[2],PlasmaZs[2]))


# plascoil_position = Observable(Point2f[])




f = Figure(); axf = Axis(f[1,1])
scpls = scatter!(plascoil_position, color=:black, markersize=20)



on(events(f).mousebutton, priority=2) do event
    global dragging, idx
    if event.button == Mouse.left
        if event.action == Mouse.press
            plt, i = pick(scpls, mouseposition_px(f), 10)
            dragging = plt == scpls
            idx = i
            @show i, dragging
            return Consume(dragging)
        elseif event.action == Mouse.release
            
            dragging = false
            return Consume(false)
        end
    end
end

on(events(f).mouseposition, priority=2) do mp
    if dragging
        plascoil_position[][idx] = mouseposition(axf)
        notify(plascoil_position)
        return Consume(true)
    end
    return Consume(false)
end



Z_vals = lift(plascoil_position) do update
    Z_coils .+ 
        [coil_fields.ψ(coil_fields.Coil(Float64(plascoil_position.val[1][1]),
            Float64(plascoil_position.val[1][2]),PlasmaArray[1].J), x,y) for x in X, y in Y]
end


contour!(axf,X,Y,Z_vals, levels=150)



f
=#





