module CoilFields

using LinearAlgebra: norm, cross, dot
using Base.Threads
using OhMyThreads: tmap
using StaticArrays
using DelimitedFiles: readdlm
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, ContinuousCallback, EnsembleProblem, EnsembleThreads, remake, terminate!
using NonlinearSolve: NonlinearProblem, solve



# include("Fourier.jl")
include("coils.jl")
include("Fields.jl")
include("io_coils.jl")
include("FieldLines.jl")
include("utilities.jl")
include("plotting.jl")



export Coil, CoilSet
export Biot_Savart, Biot_Savart!
export Biot_Savart_A, Biot_Savart_A!
export CompactLinear
export ReadCoilSet
export FieldLine
export find_axis



export plotcoil, plotcoil!
export plotcoils, plotcoils!

end # module CoilFields
