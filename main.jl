

include("./coil_fields.jl")

using GLMakie

using Contour





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



"""

"""

PlasmaRs = [3.5]
PlasmaZs = [0.0]
PlasmaJs = [-4.2e4]

PlasmaArray = [coil_fields.Coil(R,Z,J) for (R,Z,J) in zip(PlasmaRs,PlasmaZs,PlasmaJs)]


X = range(0.1,11.0,1001)
Y = range(-7.0,7.0,1001)




Z_coils = [coil_fields.ψ(CoilArray,x,y) for x in X, y in Y]

Z_plasma = [coil_fields.ψ(PlasmaArray,x,y) for x in X, y in Y]




ψ_field(r,z) = 5exp(-2(25(r-3)^2 + 5(z-0.0)^2))






# function ψ_field(r,z)
#     r_0, z_0 = (3.0,0.0)
#     a_r, a_z = (1.0,1.0)
#     α = 1
#     ω = 2

#     R = sqrt( a_r*(r-r_0-2π)^2 + a_z*(z-z_0-2π)^2 )

#     return α*tanh(ω*R) + 1
# end



function f(r,z)

    R = x
    α = [10.0,10.0]     # Slope strength
    R_0 = [-3.0, 0.0]   # Origin
    w = [0.0, 0.0]      # Width

    if R < 0
        return 1 / (1 + exp(-α * (R - R_0 - w)))
    else
        return -1 / (1 + exp(-α * (R + R_0 + w))) + 1
    end

    
end




Z_full = 0Z_coils .- [ψ_field(x,y) for x in X, y in Y]


# Z_full = Z_coils .+ Z_plasma


contour_xy = [3.5,0.0]
contour_value = coil_fields.ψ(CoilArray,contour_xy...)
contour_trace = Contour.contour(X,Y,Z_full,contour_value)
# contour_verts = hcat([[X[1],X[2]] for X in contour_trace.lines[3].vertices]...)









f = Figure(); 
axf = Axis(f[1,1])
contour!(axf,X,Y,Z_full, levels=500)
Colorbar(f[1,3],limits=(minimum(Z_full),maximum(Z_full)))
scatter!(axf,CoilRs,CoilZs,color=:red)
# scatter!(axf,PlasmaRs,PlasmaZs,color=:black)

# scatter!(axf,contour_verts[1,:], contour_verts[2,:],color=:red)

contour!(axf,X,Y,Z_full,levels=[-1e-8],color=:black)

xlims!(axf,0.1,10)
ylims!(axf,-5,5)


colsize!(f.layout, 1, Relative(0.8))
rowsize!(f.layout, 1, Relative(0.8))

axf2 = Axis(f[2,1])
lines!(axf2, X, ψ_field.(X,0.0))

axf3 = Axis(f[1,2])
lines!(axf3, Y, ψ_field.(3.0,Y))



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





