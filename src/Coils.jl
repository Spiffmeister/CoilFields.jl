module Coils

using DelimitedFiles: readdlm
using StaticArrays: SVector
using HTTP: get
using JSON: parse


include("Coils/types.jl")
include("Coils/Fourier.jl")
include("Coils/geometry.jl")
include("Coils/coils.jl")

include("Coils/read_file.jl")
include("Coils/read_quasr.jl")


export AbstractCoil, AbstractCoilGeometry, AbstractCoilSet
export Coil, CoilSet, CompositeCoilSet

export readcoilset, getcoilset

end
