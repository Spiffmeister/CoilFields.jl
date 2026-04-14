
"""
Get the `CoilSet` or `CompositeCoilSet`
"""
function GetCoilSet(location, coiltype=:quasr)
    quasrdict = get_coils_quasr(location)
    return quasr_to_coilset(quasrdict)
end



"""
Given a `quasrID` return a `JSON.Object` of the parsed dataset.
"""
function get_coils_quasr(quasrID)
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
function collect_quasr_coildata(simsopt_objs, coilkey)
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
function quasr_to_fourier(seriesval, dofs, dofnames)
    inds = findall(name -> occursin(seriesval, name), dofnames)
    series_dofs = dofs[inds]

    if occursin("c", seriesval)
        seriestype = :cos
    elseif occursin("s", seriesval)
        seriestype = :sin
    end

    return Fourier(seriestype, series_dofs, length(series_dofs))
end

"""
Given an `axis` (`"x"`, `"y"` or `"z"`) and a `Vector` of `amplitudes` and `dofnames` create a `FourierSeries`.
"""
function quasr_to_coil(axis::String, amplitudes, dofnames)
    FS = map(stype -> quasr_to_fourier(string(axis, stype), amplitudes, dofnames), ("c", "s"))
    return FourierSeries(FS...)
end

"""
 Convert a `Vector` of `amplitudes` and `dofnames` to a `FourierCurve`
"""
function quasr_to_coil(amplitudes::Vector, dofnames, current)
    x, y, z = map(ax -> quasr_to_coil(ax, amplitudes, dofnames), ["x", "y", "z"])
    FC = FourierCurve(x, y, z)
    return Coil(FC, current, length(dofnames))
end
quasr_to_coil(curvedata::NamedTuple) = quasr_to_coil(curvedata.dofs, curvedata.dofnames, curvedata.current)

"""
Take the output of `_get_coils_quasr` and convert it to a `CoilSet`
"""
function quasr_to_coilset(quasrdict)
    simsopt_objs = quasrdict["simsopt_objs"]
    # Get the coilnames and the data for each coil
    coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Coil", key)]
    # Collect all the coil data
    curve_data = [collect_quasr_coildata(simsopt_objs, coilname) for coilname in coilkeys]
    # generaate a vector of the coils
    coilset = [quasr_to_coil(curve) for curve in curve_data]
    return coilset
end
