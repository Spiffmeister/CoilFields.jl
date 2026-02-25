using Revise
using CoilFields
using HTTP
using JSON

quasrID = "0019907"

quasrdict = CoilFields.get_coils_quasr(quasrID)


function _find_quasr_coils(quasrdict)
    # Find keys with "coil" in the name
    coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Coil", key)]
    # get the fourier object which corresponds to this coil
    curvekeys = [_descend_coils(quasrdict["simsopt_objs"], key) for key in coilkeys]

    return coilkeys, curvekeys
end

function _descend_coils(simsopt_objs, coilkey)
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

    return curvekey, rotation, flip, currentkey
end

function _get_curve_data(simsopt_objs, curvekey)
    # The ID of the actual FS
    coilid = simsopt_objs[curvekey]["dofs"]["value"]

    knots = simsopt_objs[curvekey]["quadpoints"]["data"]
    series_order = simsopt_objs[curvekey]["order"]
    # The fourier amplitudes and the elements of the series
    dofs = simsopt_objs[coilid]["x"]["data"]
    dofnames = simsopt_objs[coilid]["names"]

    return coilid, knots, series_order, dofs, dofnames
end


coils_curves = _find_quasr_coils(quasrdict)

curve_data = _get_curve_data(quasrdict["simsopt_objs"], coils_curves[2][1][1])

function _quasr_to_coil(seriesval, dofs, dofnames)
    inds = findall(name -> occursin(seriesval, name), dofnames)
    series_dofs = Float64.(dofs[inds])

    if occursin("c", seriesval)
        seriestype = :cos
    elseif occursin("s", seriesval)
        seriestype = :sin
    end

    return CoilFields.Fourier(seriestype, series_dofs, length(series_dofs))
end

# seriesvals = ["x", "y", "z"], ["c", "s"]
fourierc = _quasr_to_coil("xc", curve_data[4], curve_data[5])
fouriers = _quasr_to_coil("xs", curve_data[4], curve_data[5])
x = CoilFields.FourierSeries(fourierc, fouriers)
