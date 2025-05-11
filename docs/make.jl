using Documenter, PhysicalOptics 

DocMeta.setdocmeta!(PhysicalOptics, :DocTestSetup, :(using PhysicalOptics); recursive=true)

makedocs(;
    modules=[PhysicalOptics],
    sitename="PhysicalOptics.jl",
    doctest = false,
    pages = [
        "PhysicalOptics.jl" => "index.md",
    ]
)


deploydocs(;
    repo = "github.com/JuliaPhysics/PhysicalOptics.jl.git",
    devbranch = "main",
    push_preview = true,
)
