using GLMakie
using CoilFields


coilset = ReadCoilSet("./test/coilset", skipstart=3)
coil = coilset.Coils[1]


plotcoil(coil, color=:black)

plotcoils(coilset, color=:black)
