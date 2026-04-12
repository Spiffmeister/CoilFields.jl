using Revise
using CoilFields
# using Profile

coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)
coilset = ReadCoilSet("./test/coils.W7x", :delim, endcoil_delim=iszero, endcoil_column=4, skipstart=3)

pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
B = [zeros(3) for _ in eachindex(pts)];
Biot_Savart!(B, coilset, pts, CompactLinear())
# Profile.clear_malloc_data()
# Biot_Savart!(B, coilset, pts, CompactLinear())

using BenchmarkTools
@benchmark Biot_Savart!($B, $coilset, $pts, $CompactLinear())


# quasrID = "0019907"
# quasrcoils = GetCoilSet(quasrID)

# coilfile = readdlm("./test/coilset", skipstart=3)

# using ProfileView, Cthulhu
# @profview Biot_Savart!(B, coilset, pts, CompactLinear())
# @profview Biot_Savart!(B, coilset, pts, CompactLinear())


ctypes = unique(typeof.(coilset))
if length(unique) > 1
    fme = [findall(==(c), typeof.(coilset)) for c in ctypes]
    coilgroupsets = [CoilFields.CoilGroup([coilset[Is]...]) for Is in fme]
end
find
