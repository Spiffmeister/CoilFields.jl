using Revise
using CoilFields

using GLMakie

coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)
# coilset = ReadCoilSet("./test/coils.W7x", :delim, endcoil_delim=iszero, skipstart=3, endcoil_column=4)

# Randomly generate initial points
# X₀ = [10.0, 0.0, 0.0]
# r₀ = 1.0
# pdata = CoilFields.construct_poincare(coilset, X₀, r₀, N_traj=100, t_f=1200)

# Poincare with coils
# f = plotcoils(coilset)
# scatter!(f.plot, Point3f.(pdata.points))
# f



# Choose the initial points
# X₁ = [[10.5, 0.0, 0.0], [10.9, 0.0, 0.0]]
# X₁ = [[10.5, 0.0]]
X₁ = [10.0, 0.0, 0.0]
r₀ = 0.5
pdata_exact = CoilFields.construct_poincare(coilset, X₁, r₀, N_traj=2)


# 2D poincare
scatter(Point2f.(first.(pdata_exact.points), last.(pdata_exact.points)))
# scatter!(Point2f(pdata_exact.points[1]...), color=:red)
