using Documenter, PhysicalOptics 



DocMeta.setdocmeta!(DeconvOptim, :DocTestSetup, :(using DeconvOptim); recursive=true)
makedocs(cite_bib, modules=[DeconvOptim],
         sitename="PhysicalOptics.jl",
         doctest = false,
         pages = Any[
                "PhysicalOptics.jl" => "index.md",
         ]
        )
         

deploydocs(repo = "github.com/JuliaPhysics/PhysicalOptics.jl.git")
