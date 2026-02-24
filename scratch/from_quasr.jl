using CoilFields
using HTTP
using JSON

quasrID = "0019907"

quasr_address = "https://quasr.flatironinstitute.org/simsopt_serials/"
quasrJSON = HTTP.get(
    string(quasr_address, quasrID[1:4], "/", "serial", quasrID, ".json")
)
quasrstr = String(quasrJSON.body)

quasrdict = JSON.parse(quasrstr, allownan=true)

coilkeys = [key for key in keys(quasrdict["simsopt_objs"]) if occursin("Curve", key)]

# Select a coil object
coilobj = quasrdict["simsopt_objs"][coilkeys[1]]
knots = coilobj["quadpoints"]["data"]
order = coilobj["order"]
dofs_key = coilobj["dofs"]["value"]

# Get the dofs of the coil object and the dof names
dofs = quasrdict["simsopt_objs"][dofs_key]["x"]["data"]
dofnames = quasrdict["simsopt_objs"][dofs_key]["names"]

# How the series are generally called
seriesvals = ["x", "y", "z"], ["c", "s"]
for ax in seriesvals[1]
    for stype in seriesvals[2]
        inds = findall(name -> occursin(string(ax, stype), name), dofnames)
        series_dofs = Float64.(dofs[inds])

        if stype == "c"
            seriestype = :cos
        elseif stype == "s"
            seriestype = :sin
        end


    end
end
