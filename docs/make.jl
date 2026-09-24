using Documenter
using CubatureRules

# Julia's Markdown splits a table row on every `|`, even inside code, so `e^{-|x|²}` in a
# cell silently adds columns and the table renders broken without any warning. Refuse to
# build instead.
let bad = String[]
    for (root, _, files) in walkdir(joinpath(@__DIR__, "src")), f in files
        endswith(f, ".md") || continue
        lines = readlines(joinpath(root, f))
        i = 1
        while i <= length(lines)
            j = i
            while j <= length(lines) && startswith(lines[j], "|")
                j += 1
            end
            n = [count(==('|'), l) for l in lines[i:j-1]]
            isempty(n) || allequal(n) || push!(bad, "$(relpath(joinpath(root, f), @__DIR__)):$i")
            i = max(j, i + 1)
        end
    end
    isempty(bad) || error("tables whose rows have different numbers of cells (a `|` inside a cell?): " *
                          join(bad, ", "))
end

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
            "I want an error estimate" => "tutorial/errors.md",
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
            "Error estimates" => "design/errors.md",
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
