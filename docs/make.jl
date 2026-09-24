using Documenter
using CubatureRules

makedocs(;
    sitename = "CubatureRules.jl",
    modules = [CubatureRules],
    # the API reference is one page on purpose, so that it can be searched with Ctrl-F
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                             size_threshold_warn = 150 * 2^10),
    # every public name must be in the manual; internal helpers keep their docstrings for
    # readers of the source without being listed
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        "Tutorial" => [
            "I just want a rule" => "tutorial/first-rule.md",
            "I want to choose a rule" => "tutorial/choosing.md",
            "I want arbitrary precision" => "tutorial/precision.md",
            "I want to integrate to a tolerance" => "tutorial/tolerance.md",
            "I want my own domain or mesh" => "tutorial/domains.md",
            "I want an unusual weight" => "tutorial/weights.md",
            "I want to verify a rule" => "tutorial/verifying.md",
            "I want to cite a rule" => "tutorial/citing.md",
            "I want performance" => "tutorial/performance.md",
            "I want to build a family" => "tutorial/families.md",
            "I want to provide external data" => "tutorial/providers.md",
            "Troubleshooting" => "tutorial/troubleshooting.md",
        ],
        "Architecture & Design" => [
            "The pipeline" => "design/pipeline.md",
            "Domains and measures" => "design/domains.md",
            "Precision and guard digits" => "design/precision.md",
            "Exactness claims" => "design/claims.md",
            "Certificates" => "design/certificates.md",
            "Verification" => "design/verification.md",
            "Selection" => "design/selection.md",
            "Seed strategies" => "design/seeds.md",
            "Orbit algebra" => "design/symmetry.md",
            "Measures given by moments" => "design/moments.md",
            "Application and transport" => "design/application.md",
            "Sequences and adaptivity" => "design/adaptive.md",
            "Provenance and licensing" => "design/provenance.md",
            "External providers" => "design/providers.md",
        ],
        "Catalogue" => [
            "Overview" => "catalogue/index.md",
            "Interval" => "catalogue/interval.md",
            "Simplex" => "catalogue/simplex.md",
            "Box" => "catalogue/box.md",
            "Sphere" => "catalogue/sphere.md",
            "Ball" => "catalogue/ball.md",
            "Unbounded" => "catalogue/unbounded.md",
        ],
        "API" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(; repo = "github.com/IlianPihlajamaa/CubatureRules.jl.git", push_preview = false)
