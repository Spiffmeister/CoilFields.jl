using CoilFields

using FixedSizeArrays

base_coilset = ReadCoilSet("./test/coilset", :delim, skipstart=3)


ccg = [[FixedSizeArray{eltype(base_coilset)}(undef, 3) for _ in 1:base_coilset[1].length] for _ in 1:length(base_coilset.Coils)]
for i_coil in 1:length(base_coilset.Coils)
    for j_pt in 1:length(base_coilset.Coils[i_coil].Geometry)
        ccg[i_coil][j_pt] .= base_coilset.Coils[i_coil].Geometry[j_pt]
    end
end

coils = [Coil(ccg[i_coil], base_coilset[i_coil].J, base_coilset[i_coil].length) for i_coil in 1:length(base_coilset.Coils)]
