using Revise
using CoilFields
using BenchmarkTools

coilset = ReadCoilSet("./test/coils.W7x", :delim, endcoil_delim=iszero, endcoil_column=4, skipstart=3)
pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
@benchmark Biot_Savart($coilset, $pts, $CompactLinear())


# using ProfileView
# using Cthulhu

# @profview Biot_Savart(coilset, pts, CompactLinear())
# @profview Biot_Savart(coilset, pts, CompactLinear())
