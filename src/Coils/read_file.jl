"""
    readcoilset(filename, coiltype=:delim; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"], endcoil_column=0)

For reading in delimited coil files.



Optional Inputs:
- endcoil_delim: What string is used to mark the final row for a coil
    Can also be a `Function` (i.e. `iszero`)
- skipstart: the number of rows to skip in the read in, default `0`
- filelayout: the column labels in the file
- endcoil_column: The column label which corresponds to the end of a coil.


As an example, consider reading the coilset included in the `test` folder:
```julia
coilset = readcoilset("./test/coilset", skipstart=3)
```

Returns a [`CoilSet`](@ref) if all `Coil` objects are of the same type or [`CompositeCoilSet`](@ref) otherwise.
"""
function readcoilset(filename, coiltype=:delim; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"], endcoil_column=0)

    if coiltype == :delim
        return read_coils_mod(filename, skipstart, filelayout, endcoil_delim, endcoil_column)
    end
end



"""
Read coils from a file with the coil delimiter `mod`
"""
function read_coils_mod(filename, skipstart, filelayout, endcoil_delim, endcoil_column=0)

    fread = readdlm(filename, skipstart=skipstart)

    idx_startofcoil = 1
    endcoil_column < 1 ? endcoil_column = length(fread[1, :]) : nothing
    if typeof(endcoil_delim) <: String
        nextcoil(str) = occursin(endcoil_delim, str)
    elseif typeof(endcoil_delim) <: Function
        nextcoil = endcoil_delim
    end
    ncoils = sum(nextcoil.(fread[:, endcoil_column]))
    coils = []
    # coiltype = []

    for _ in 1:ncoils
        # Scan though the lines and find the
        @views idx_endofcoil = findnext(nextcoil, fread[:, endcoil_column], idx_startofcoil + 1)

        coil = Coil(
            Tuple(SVector{3}(Float64.(row[1:3])) for row in eachrow(fread[idx_startofcoil:idx_endofcoil, :])),
            fread[idx_startofcoil, 4],
            idx_endofcoil - idx_startofcoil + 1
        )
        push!(coils, coil)
        idx_startofcoil = idx_endofcoil + 1
    end


    # Get the coil groups by finding unique coil types
    ctypes = unique(typeof.(coils))
    cinds = [findall(==(c), typeof.(coils)) for c in ctypes]
    coilgroupsets = [CoilSet([coils[Is]...]) for Is in cinds]
    if length(coilgroupsets) == 1
        coilset = coilgroupsets[1]
    else
        coilset = CompositeCoilSet(Tuple(coilgroupsets))
    end

    return coilset

end
