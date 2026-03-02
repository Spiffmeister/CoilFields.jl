using Revise
using CoilFields
using HTTP
using JSON

quasrID = "0019907"

# quasrdict = CoilFields._get_coils_quasr(quasrID)




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

function _quasr_to_coil(amplitudes::Vector, labels, current)
    x, y, z = map(ax -> _quasr_to_coil(ax, amplitudes, labels), ["x", "y", "z"])
    FC = CoilFields.FourierCurve(x, y, z)
    return CoilFields.Coil(FC, current, length(labels))
end

_quasr_to_coil(curvedata::NamedTuple) = _quasr_to_coil(curvedata.dofs, curvedata.dofnames, curvedata.current)


function _quasr_to_coilset(quasrdict)
    simsopt_objs = quasrdict["simsopt_objs"]
    # Get the coilnames and the data for each coil
    coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Coil", key)]
    # Collect all the coil data
    curve_data = [_collect_quasr_coildata(simsopt_objs, coilname) for coilname in coilkeys]
    # generaate a vector of the coils
    coilset = [_quasr_to_coil(curve) for curve in curve_data]
    return CoilFields.CoilSet(coilset)
end


# curve_data = _get_curve_data(quasrdict["simsopt_objs"], coils_curves[2][1][1])
# seriesvals = ["x", "y", "z"], ["c", "s"]
# x = _quasr_to_coil("x", curve_data[4], curve_data[5])
# f = _quasr_to_coil(curve_data.dofs, curve_data.dofnames)
coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Coil", key)]
curve_data = [_collect_quasr_coildata(quasrdict["simsopt_objs"], coilname) for coilname in coilkeys]


ooga = _quasr_to_coilset(quasrdict)
