"""
    ReadCoilSet(filename; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"])

For reading in the most common format of coil file.

Optional Inputs:
- endcoil_delim: What string is used to mark the final row for a coil
- skipstart: the number of rows to skip in the read in
- filelayout: the column labels in the file


Endcoil delim:
    - Supports `(function, colindex)`

Returns a [`CoilSet`](@ref) object.
"""
function ReadCoilSet(filename, coiltype; endcoil_delim="mod", skipstart=0, filelayout=["x", "y", "z", "current"], endcoil_column=0)

    if coiltype == :delim
        return _read_coils_mod(filename, skipstart, filelayout, endcoil_delim, endcoil_column)
    end
end


function GetCoilSet(location, coiltype=:quasr)
    quasrdict = _get_coils_quasr(location)
    return _quasr_to_coilset(quasrdict)
end







"""
Read coils from a file with the coil delimiter `mod`
"""
function _read_coils_mod(filename, skipstart, filelayout, endcoil_delim, endcoil_column=0)

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
    coilgroupsets = [CoilFields.CoilSet([coils[Is]...]) for Is in cinds]
    if length(coilgroupsets) == 1
        coilset = coilgroupsets[1]
    else
        coilset = CompositeCoilSet(Tuple(coilgroupsets))
    end

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


"""
Gather all the information needed to construct a coil from the QUASR database

`simopt_objs` is the `simsopt_dict['simsopt_objs']` dictionary and `coilkey` is a coil ID in the `simsopt_objs`
"""
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

    # Need to descend coils until we find the correct object
    function getkey(currentname, depth=10)
        its = 1
        while occursin("ScaledCurrent", currentname) & its < depth
            currentname = simsopt_objs[currentname]["current_to_scale"]["value"]
            if !occursin("ScaledCurrent", currentname)
                return currentname
            end
            its += 1
        end
    end
    currentname = getkey(currentkey)
    currentid = simsopt_objs[currentname]["dofs"]["value"]
    current = Float64.(simsopt_objs[currentid]["x"]["data"])[1]
    # currentscale = Float64.(simsopt_objs[currentkey])

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
        dofnames=dofnames,
        current=current
    )

    return coildata
end


"""
Convert a `Vector` of `dofs` to a `Fourier` series from the QUASR database.

`seriesval` is used to search the `dofnames` for the given series to create.
"""
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

"""
Given an `axis` (`"x"`, `"y"` or `"z"`) and a `Vector` of `amplitudes` and `dofnames` create a `FourierSeries`.
"""
function _quasr_to_coil(axis::String, amplitudes, dofnames)
    FS = map(stype -> _quasr_to_fourier(string(axis, stype), amplitudes, dofnames), ("c", "s"))
    return CoilFields.FourierSeries(FS...)
end

"""
 Convert a `Vector` of `amplitudes` and `dofnames` to a `FourierCurve`
"""
function _quasr_to_coil(amplitudes::Vector, dofnames, current)
    x, y, z = map(ax -> _quasr_to_coil(ax, amplitudes, dofnames), ["x", "y", "z"])
    FC = CoilFields.FourierCurve(x, y, z)
    return CoilFields.Coil(FC, current, length(dofnames))
end
_quasr_to_coil(curvedata::NamedTuple) = _quasr_to_coil(curvedata.dofs, curvedata.dofnames, curvedata.current)

"""
Take the output of `_get_coils_quasr` and convert it to a `CoilSet`
"""
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
