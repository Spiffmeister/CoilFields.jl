# # Field line tracing

using CoilFields, GLMakie;
coilset = ReadCoilSet("../../../test/coilset", skipstart=3);

# Poincare section create currently only works on a single slice, multiple slices is too be implemented.
#
# Randomly generate initial points within an ``r_0`` ball of ``X_0``
X₀ = [10.0, 0.0, 0.0]
r₀ = 1.0
pdata = construct_poincare(coilset, X₀, r₀, N_traj=100, event=true)
scatter(Point2f.(first.(pdata.points), last.(pdata.points)))


# Otherwise we can choose the initial points  explicitly
X₁ = [[10.5, 0.0, 0.0], [10.9, 0.0, 0.0]]
pdata = construct_poincare(coilset, X₁, event=true)


# 2D poincare
scatter(Point2f.(first.(pdata.points), last.(pdata.points)))
