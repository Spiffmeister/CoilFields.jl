module CoilFields

using LinearAlgebra: norm, cross, dot
using Base.Threads
using Folds

include("Fourier.jl")
include("coils.jl")


export Coil, CoilSet
export Biot_Savart, Biot_Savart!
export CompactLinear

end # module CoilFields
