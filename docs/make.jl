
push!(LOAD_PATH,"../src/")

using Documenter, Literate, CoilFields

LitPathExample = joinpath(@__DIR__,"..","examples","reading.jl")
DocSrc = joinpath(@__DIR__,"src","examples") #.md creation path
Literate.markdown(LitPathExample,DocSrc,codefence="```text" => "```")


makedocs(sitename="CoilFields",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "examples/reading.md"
        ]
    ],
    modules=[CoilFields],
    format=Documenter.HTML(prettyurls=false),
    warnonly = Documenter.except(:linkcheck,:footnote)
    )

# deploydocs(
#     repo = "github.com/Spiffmeister/CoilFields.jl.git",
# )
