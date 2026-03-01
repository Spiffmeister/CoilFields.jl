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





"""
Given a `quasrID` return a `JSON.Object` of the parsed dataset.
"""
function _get_coils_quasr(quasrID)
    quasr_address = "https://quasr.flatironinstitute.org/simsopt_serials/"
    quasrJSON = get(
        string(quasr_address, quasrID[1:4], "/", "serial", quasrID, ".json")
    )
    quasrstr = String(quasrJSON.body)

    quasrdict = parse(quasrstr, allownan=true)
    return quasrdict
end


function _collect_quasr_coildata(simsopt_objs, coilkey)
    # A curve is either a
    curvekey = simsopt_objs[coilkey]["curve"]["value"]
    currentkey = simsopt_objs[coilkey]["current"]["value"]
    rotation = zero(Float64)
    flip = false
    # Its possible we don't get the fourier object
    if occursin("RotatedCurve", curvekey)
        rotation = Float64(simsopt_objs[curvekey]["phi"])
        flip = simsopt_objs[curvekey]["flip"]
        # we can now replace curvekey that we have what we need
        curvekey = simsopt_objs[curvekey]["curve"]["value"]
    end

    coilid = simsopt_objs[curvekey]["dofs"]["value"]
    knots = simsopt_objs[curvekey]["quadpoints"]["data"]
    order = Int(simsopt_objs[curvekey]["order"])

    dofs = Vector{Float64}(simsopt_objs[coilid]["x"]["data"])
    dofnames = Vector{String}(simsopt_objs[coilid]["names"])

    coildata = (curvename=curvekey,
        rotation=rotation,
        flipped=flip,
        currentname=currentkey,
        coilid=coilid,
        knots=knots,
        order=order,
        dofs=dofs,
        dofnames=dofnames
    )

    return coildata
end


function _quasr_to_fourier(seriesval, dofs, dofnames)
    inds = findall(name -> occursin(seriesval, name), dofnames)
    series_dofs = dofs[inds]

    if occursin("c", seriesval)
        seriestype = :cos
    elseif occursin("s", seriesval)
        seriestype = :sin
    end

    return CoilFields.Fourier(seriestype, series_dofs, length(series_dofs))
end

function _quasr_to_coil(axis::String, amplitudes, labels)
    FS = map(stype -> _quasr_to_fourier(string(axis, stype), amplitudes, labels), ("c", "s"))
    return CoilFields.FourierSeries(FS...)
end

function _quasr_to_coil(amplitudes::Vector, labels)
    x, y, z = map(ax -> _quasr_to_coil(ax, amplitudes, labels), ["x", "y", "z"])
    return CoilFields.FourierCurve(x, y, z)
end

_quasr_to_coil(curvedata::NamedTuple) = _quasr_to_coil(curvedata.dofs, curvedata.dofnames)


function _quasr_to_coilset(quasrdict)
    simsopt_objs = quasrdict["simsopt_objs"]
    # Get the coilnames and the data for each coil
    coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Coil", key)]
    # Collect all the coil data
    curve_data = [_collect_quasr_coildata(simsopt_objs, coilname) for coilname in coilkeys]
    # generaate a vector of the coils
    coilset = [_quasr_to_coil(curve) for curve in curve_data]
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
