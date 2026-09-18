# Regenerate test/reference_hashes.toml. Only run this when a rule is changed on purpose,
# and say why in the commit message: the file is the reproducibility contract.
using CubatureRules, TOML
include(joinpath(@__DIR__, "reference_cases.jl"))

hashes = Dict(name => rule_hash(make()) for (name, make) in REFERENCE_CASES)
open(joinpath(@__DIR__, "reference_hashes.toml"), "w") do io
    println(io, "# SHA-256 content hashes of reference rules; see test/test_reproducibility.jl.")
    TOML.print(io, hashes; sorted = true)
end
