module plotting

using Makie
using CoilFields

import CoilFields: plotcoil, plotcoil!
import CoilFields: plotcoils, plotcoils!


Makie.@recipe PlotCoil begin
    color = :black
end

Makie.preferred_axis_type(::PlotCoil) = LScene



function Makie.plot!(input::PlotCoil{<:Tuple{<:CoilFields.Coil}})
    coil_pts = Makie.Point3f.(input[1][].Geometry)
    Makie.lines!(input, [coil_pts...], color=input.color)
    input
end



Makie.@recipe PlotCoils begin
    color = :black
end

Makie.preferred_axis_type(::PlotCoils) = LScene

function Makie.plot!(input::PlotCoils{<:Tuple{<:CoilFields.CoilSet}})
    for coil in input[1][].Coils
        coil_pts = Makie.Point3f.(coil.Geometry)
        Makie.lines!(input, [coil_pts...], color=input.color)
    end
    input
end


end
