using Revise

using CoilFields

include("Frames.jl")
using .Frames

m_max = 3
n_max = 3

X1 = zeros(m_max, n_max)
X1[1, 1] = 3.0
X1[2, 1] = 1.0
X1[2, 2] = 0.4

X1_amn = (
    (3.0, 0, 0),
    (1.0, 1, 0),
    (0.4, 1, 1)
)

X2_amn = (
    (1.0, 1, 0),
    (-0.4, 1, 1),
    (-0.25, 0, 1)
)

X¹ = CoilFields.Coils.Fourier2D(:cos, X1_amn, 3)
X² = CoilFields.Coils.Fourier2D(:sin, X2_amn, 3)





## # Initialise a circular GFrame as in workshop materials
# The origin curve
X₀ = CoilFields.Coils.FourierCurve(
    x=CoilFields.Coils.FourierSeries(
        CoilFields.Coils.Fourier(:cos, 5.0, 1),
        CoilFields.Coils.Fourier(:sin, 5.0, -1)
    ))

# N points outwards
N = CoilFields.Coils.FourierCurve(
    x=CoilFields.Coils.FourierSeries(
        CoilFields.Coils.Fourier(:cos, 1.0, 1),
        CoilFields.Coils.Fourier(:sin, 1.0, -1)
    ))

# B is a vector pointing in the z direction
B(q) = (0.0, 0.0, 1.0)


X₀(0.0)


h = Frames.GFrame(X₀, N, B)


q = zeros(3)

Frames.apply_frame(q, h)
