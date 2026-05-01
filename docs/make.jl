
push!(LOAD_PATH, "../src/")

using Documenter, Literate, CoilFields

LitPathExampleFLT = joinpath(@__DIR__, "..", "examples", "fieldlinetracing.jl")
LitPathExamplePlot = joinpath(@__DIR__, "..", "examples", "plots.jl")
DocSrc = joinpath(@__DIR__, "src", "examples") #.md creation path
Literate.markdown(LitPathExamplePlot, DocSrc)
Literate.markdown(LitPathExampleFLT, DocSrc)



makedocs(sitename="CoilFields",
    pages=[
        "Home" => "index.md",
        "Reading Coils" => [
            "reading.md"
        ],
        "Coil sets" => [
            "coils.md",
            "examples/plots.md"
        ],
        "Magnetic field" => [
            "magneticfield.md",
            "examples/fieldlinetracing.md"
        ]
    ],
    modules=[CoilFields],
    format=Documenter.HTML(prettyurls=false),
    warnonly=Documenter.except(:linkcheck, :footnote)
)

deploydocs(
    repo = "github.com/Spiffmeister/CoilFields.jl.git",
)
