using Revise
using CoilFields
# using Profile

coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)

pts = [[0.0, 0.0, z] for z in range(-1.0, 1.0, 1_000)];
B = [zeros(3) for _ in eachindex(pts)];
Biot_Savart!(B, coilset, pts, CompactLinear())

# quasrID = "0019907"
# quasrcoils = GetCoilSet(quasrID)
