using CoilFields
using GLMakie

# 
# We provide some utilities for plotting coil sets in 3D
# 
coilset = ReadCoilSet("./test/coilset", skipstart=3)

#
# A single coil object can be plotted with
coil = coilset[1].Coils[1]
plotcoil(coil, color=:black)


# A coilset can be plotted with
plotcoils(coilset, color=:black)
