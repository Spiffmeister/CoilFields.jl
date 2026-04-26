module CoilFields

using LinearAlgebra: norm, cross, dot
using OhMyThreads: tmap
using StaticArrays
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, ContinuousCallback, EnsembleProblem, EnsembleThreads, remake, terminate!
using NonlinearSolve: NonlinearProblem, solve


include("constants.jl")

# include("coils.jl")
include("Coils.jl")
using .Coils

include("Fields.jl")

# include("io_coils.jl")
include("FieldLines.jl")
include("utilities.jl")
include("plotting.jl")



export Coil, CoilSet, CompositeCoilSet
export biot_savart, biot_savart!
# export biot_savart_A, biot_savart_A!
export CompactLinear
export FieldLine
export find_axis

export readcoilset
export getcoilset

export construct_poincare


export plotcoil, plotcoil!
export plotcoils, plotcoils!

end # module CoilFields
