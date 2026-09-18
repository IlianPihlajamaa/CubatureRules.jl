# Every file in src/data must have an entry in src/data/PROVENANCE.toml recording its
# author, source, licence and retrieval date (PLAN §0.2). This test is the CI check.
using Test, TOML

datadir = joinpath(@__DIR__, "..", "src", "data")
manifest = TOML.parsefile(joinpath(datadir, "PROVENANCE.toml"))
entries = get(manifest, "file", Any[])
recorded = Set(e["path"] for e in entries)

for f in readdir(datadir)
    f == "PROVENANCE.toml" && continue
    @test f in recorded
end
required = ["path", "description", "authors", "source", "license", "retrieved", "citation"]
for e in entries
    for k in required
        @test haskey(e, k) && !isempty(string(e[k]))
    end
    @test isfile(joinpath(datadir, e["path"]))
    @test !occursin("quadpy", lowercase(string(e["source"])))   # never transcribed from quadpy
end
