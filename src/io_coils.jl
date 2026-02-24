"""
    ReadCoilSet(filename; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"])

For reading in the most common format of coil file.

Optional Inputs:
- endcoil_delim: What string is used to mark the final row for a coil
- skipstart: the number of rows to skip in the read in
- filelayout: the column labels in the file

Returns a [`CoilSet`](@ref) object.
"""
function ReadCoilSet(filename, coiltype; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"])

    if coiltype == :delim
        return _read_coils_mod(filename, skipstart, filelayout)
    end
end


"""
Read coils from a file with the coil delimiter `mod`
"""
function _read_coils_mod(filename, skipstart, filelayout)

    endcoil_delim = "mod"

    fread = readdlm(filename, skipstart=skipstart)

    ncoils = sum(occursin.(endcoil_delim, fread[:, end]))

    idx_startofcoil = 1

    coils = []

    for i_coil in 1:ncoils
        idx_endofcoil = findnext(occursin.(endcoil_delim, fread[:, end]), idx_startofcoil + 1)
        coil = Coil(
            Tuple(SVector{3}(Float64.(row[1:3])) for row in eachrow(fread[idx_startofcoil:idx_endofcoil, :])),
            fread[idx_startofcoil, 4],
            idx_endofcoil - idx_startofcoil + 1
        )
        push!(coils, coil)
        idx_startofcoil = idx_endofcoil + 1
    end

    coilset = CoilSet([coils...])

    return coilset
end




function WriteCoilSet(filename, coilset; format=:hdf5, header=nothing)

    if format == :ascii
        _write_coil_ascii(filename, coilset, header)
    elseif format == :hdf5
    else
        error("format must be :hdf5 or :ascii")
    end

end

function _write_coil_hdf5(filename, coilset)
    error("Not implemented")
end

function _write_coil_ascii(filename, coilset, header=nothing)

    open(filename, "w") do io
        if !isnothing(header)
            for header_row in eachrow(header)
                writedlm(io, header_row, "")
            end
        end

        # iterate over the coils. coilnumber is required for knowing when a new coil starts in the "standard" file format
        for (coilnumber, coil) in enumerate(coilset)
            for (i, xyz) in enumerate(coil.Geometry)
                # The standard file format has J=0 at the last index and a 'COILNUMBER mod_PADDEDCOILNUMBER'
                #   where PADDEDCOILNUMBER is just the COILNUMBER padded to the left with zeros so that it is
                #   three characters long
                if i ≠ lastindex(coil.Geometry)
                    writedlm(io, [xyz, coil.J], '\t')
                else
                    writedlm(io, [xyz, 0.0, "coilnumber mod_$(lpad(coilnumber,3,"0"))"], '\t')
                end
            end
        end
    end

end
