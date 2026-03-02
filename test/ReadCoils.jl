using CoilFields

using DelimitedFiles

coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)


quasrID = "0019907"
quasrcoils = GetCoilSet(quasrID)

# coilfile = readdlm("./test/coilset", skipstart=3)
