using CoilFields
using LinearAlgebra


# Test from coilset
coilset = readcoilset("./test/coilset", :delim, skipstart=3)

X = [10.0, 0.0]

axis_origin = find_axis(X, coilset)

# axis_expected = [10.00949838182673, 0.0, 8.133575801395752e-7]
axis_expected = [10.00949838182673, 0.0, 2.9179442361179924e-7]

# Adjust tolerances/method for more accuracy
@test norm(axis_expected - axis_origin) ≤ 1e-6
