using SafeTestsets



@safetestset "Biot savart" begin
    include("BiotSavart.jl")
end

@safetestset "Reading test" begin
    include("ReadCoils.jl")
end

@safetestset "Axis finding" begin
    include("axis_finder.jl")
end
