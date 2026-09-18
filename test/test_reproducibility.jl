# Bitwise reproducibility (PLAN §5 item 9): identical `rule(...)` calls must produce
# identical output across sessions, Julia versions and platforms. The reference hashes in
# reference_hashes.toml are checked on every CI job, across the whole OS × version matrix.
#
# Regenerate (only when a rule is intentionally changed) with
#     julia --project test/reference_hashes.jl
using CubatureRules, Test, TOML

include("reference_cases.jl")

ref = TOML.parsefile(joinpath(@__DIR__, "reference_hashes.toml"))
for (name, make) in REFERENCE_CASES
    r = make()
    # repeat construction in the same session: identical
    @test rule_hash(r) == rule_hash(make())
    @test haskey(ref, name)
    if haskey(ref, name)
        h = rule_hash(r)
        h == ref[name] || @info "hash mismatch" name h expected = ref[name]
        @test h == ref[name]
    end
end
