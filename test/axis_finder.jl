using CoilFields



# Test from coilset
coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)

X = [10.0, 0.0]

axis_origin = find_axis(X, coilset)

axis_expected = [10.00949838182673, 0.0, 8.133575801395752e-7]
