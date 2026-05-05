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
# X² = CoilFields.Coils.Fourier2D(:sin, X2_amn, 3)

X¹(0.0, 0.0)




Xa = [1.0, 1.0]
Xm = 1

X = CoilFields.Coils.Fourier(:cos, Xa, Xm)
N = CoilFields.Coils.Fourier(:cos, Xa, Xm)
B = CoilFields.Coils.Fourier(:cos, Xa, Xm)






h = Frames.GFrame(X, N, B)


# q =

# apply_frame()
