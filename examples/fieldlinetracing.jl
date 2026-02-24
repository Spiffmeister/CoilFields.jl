using Revise
using CoilFields

using GLMakie

coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)

# Randomly generate initial points
X₀ = [10.0, 0.0, 0.0]
r₀ = 1.0
pdata = CoilFields.construct_poincare(coilset, X₀, r₀, N_traj=100, t_f=1200)

# Poincare with coils
# f = plotcoils(coilset)
# scatter!(f.plot, Point3f.(pdata.points))
# f



# Choose the initial points
X₁ = [[10.5, 0.0, 0.0], [10.9, 0.0, 0.0]]
pdata_exact = CoilFields.construct_poincare(coilset, X₁, nothing, t_f=2400)


# 2D poincare
# scatter(Point2f.(first.(pdata_exact.points), last.(pdata_exact.points)))
