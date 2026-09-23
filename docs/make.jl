using Documenter
using CubatureRules

makedocs(;
    sitename = "CubatureRules.jl",
    modules = [CubatureRules],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "Design" => "design.md",
        "Rule families" => "families.md",
        "Weights given by moments" => "moments.md",
        "Delegated families" => "delegation.md",
        "Verification" => "verification.md",
        "Adding a family" => "extending.md",
        "API" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(; repo = "github.com/IlianPihlajamaa/CubatureRules.jl.git", push_preview = false)
