# # Plotting coils

using CoilFields, GLMakie;
coilset = ReadCoilSet("../../../test/coilset", skipstart=3);

# There are some utilities for plotting coil sets in 3D using any [`Makie`](https://docs.makie.org/stable/) extension (for instance GLMakie).
# A single coil object can be plotted with

coil = coilset[1]
plotcoil(coil, color=:black)


# A coilset can be plotted with

plotcoils(coilset, color=:black)
