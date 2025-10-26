module CoilFields

using LinearAlgebra: norm, cross, dot
using Base.Threads
using OhMyThreads: tmap
using DelimitedFiles: readdlm

# include("Fourier.jl")
include("coils.jl")
include("ReadCoils.jl")
include("FieldLines.jl")
include("plotting.jl")

export Coil, CoilSet
export Biot_Savart, Biot_Savart!
export Biot_Savart_A, Biot_Savart_A!
export CompactLinear
export ReadCoilSet
export FieldLine

export plotcoil, plotcoil!
export plotcoils, plotcoils!

end # module CoilFields
