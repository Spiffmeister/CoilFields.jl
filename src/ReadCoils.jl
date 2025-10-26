function ReadCoilSet(filename; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"])

    fread = readdlm(filename, skipstart=skipstart)

    ncoils = sum(occursin.(endcoil_delim, fread[:, end]))

    idx_startofcoil = 1

    coils = []

    for i_coil in 1:ncoils
        idx_endofcoil = findnext(occursin.(endcoil_delim, fread[:, end]), idx_startofcoil + 1)
        coil = Coil(
            [Float64.(row[1:3]) for row in eachrow(fread[idx_startofcoil:idx_endofcoil, :])],
            fread[idx_startofcoil, 4],
            idx_endofcoil - idx_startofcoil + 1
        )
        push!(coils, coil)
        idx_startofcoil = idx_endofcoil + 1
    end

    coilset = CoilSet([coils...])

    return coilset
end
