# Interop with Lebedev.jl, which is GPL-3. Loading it here is *use*, not distribution: no
# part of it is redistributed by this package, and the test environment is not shipped.
#
# Plain `using` on both: the family here is `LebedevRule`, so the package name `Lebedev` is
# free. That the two coexist is part of what this file checks.
using CubatureRules, Test
import CubatureRules: selectable   # public, not exported
using Lebedev
const CR = CubatureRules

@testset "one family, two sources of seeds" begin
    @test UpstreamLebedev === CR.LebedevRule{CR.LebedevJLSeeds}
    @test LebedevRule() isa LebedevRule            # the in-house default
    @test LebedevRule() === CR.LebedevRule{CR.InHouseSeeds}()
    @test CR.family_name(UpstreamLebedev()) == CR.family_name(LebedevRule()) == "Lebedev"
    # they are told apart by the description, which carries the licence
    @test CR.describe_family(LebedevRule()) == "Lebedev"
    @test occursin("GPL-3", CR.describe_family(UpstreamLebedev()))
    @test CR.family_license(LebedevRule()) == ""
    @test occursin("GPL-3.0", CR.family_license(UpstreamLebedev()))
    @test CR.missing_dependency(UpstreamLebedev()) === nothing      # the extension is loaded
end

@testset "available lists both; rule takes only one" begin
    names = [row.family for row in available(Sphere{3}(); degree = 11)]
    @test "Lebedev" in names                       # ours
    @test "Lebedev (Lebedev.jl, GPL-3)" in names   # theirs, with the terms in the name
    # the selector will not take theirs on its own
    @test CR.selectable(LebedevRule()) === true
    @test CR.selectable(UpstreamLebedev()) === false
    @test all(c -> selectable(c.family), CR.gather(Sphere{3}(), 19, Float64))
    @test any(c -> !selectable(c.family), CR.gather(Sphere{3}(), 19, Float64; all = true))
    # at degree 19 our table stops, so the default is the product rule — with a warning,
    # because a smaller rule existed and was passed over for its licence
    r = @test_logs (:warn, r"licence is not one") match_mode=:any rule(Sphere{3}(); degree = 19)
    @test family(r) == "SphereProduct"
    # asked for explicitly, the smaller rule comes back
    r2 = rule(Sphere{3}(); degree = 19, copyleft = true)
    @test family(r2) == "Lebedev"
    @test npoints(r2) == 146
    @test npoints(r2) < npoints(r)
    @test occursin("copyleft = true", provenance(r2).selection)
    # the warning can be switched off
    CR.license_warnings!(false)
    @test_logs rule(Sphere{3}(); degree = 19)          # no warning at all
    CR.license_warnings!(true)
end

@testset "upstream rules carry their terms" begin
    for (d, n) in ((3, 6), (11, 50), (29, 302))
        r = rule(UpstreamLebedev(), Sphere{3}(); degree = d)
        v = check(r)
        @test npoints(r) == n
        @test sum(weights(r)) ≈ 4π rtol = 1e-13        # rescaled: Lebedev.jl normalises to 1
        @test v.exact && v.sharp === true && v.positive && v.symmetric === true
        @test occursin("GPL-3.0", provenance(r).license)
        @test occursin("not covered by CubatureRules.jl's MIT licence", provenance(r).license)
        @test occursin("Lebedev.jl", provenance(r).seed_source)
        @test occursin("LebedevLaikov1999", cite(r))   # the citation they ask for
    end
    # Float64 is a pass-through: nothing is solved, so nothing is certified
    r = rule(UpstreamLebedev(), Sphere{3}(); degree = 11)
    @test certificate(r) === nothing
    @test occursin("delivered as it stands", join(provenance(r).path, " "))
    @test npoints(UpstreamLebedev(), Sphere{3}(), 125) == 5294
    @test_throws NoRuleError rule(UpstreamLebedev(), Sphere{3}(); degree = 131)
end

@testset "refining an upstream rule keeps its licence" begin
    # refining is allowed: a caller may work with data they hold. What travels with the
    # result is the licence, because the refined rule is a derived work of their table.
    r = rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 40)
    @test eltype(r) === BigFloat
    @test npoints(r) == 302
    v = check(r)
    @test v.exact && v.sharp === true && v.positive && v.symmetric === true
    dev = CR.with_bits(() -> abs(sum(weights(r)) - 4 * BigFloat(π)), CR.digits_to_bits(40) + 64)
    @test dev < big(10.0)^-38
    @test occursin("GPL-3.0", provenance(r).license)             # unchanged by refinement
    @test occursin("derived work", join(provenance(r).path, " "))
    @test occursin("structure recovered", join(provenance(r).path, " "))
    # and it is certified like any refined rule
    @test certificate(r) !== nothing
    @test certificate(r).residual < big(10.0)^-38
    @test occursin("invariants", certificate(r).equations)
end

@testset "recovering an orbit structure from bare points" begin
    # what makes refining someone else's table possible at all
    for ord in (11, 29, 47)
        x, y, z, w = lebedev_by_order(ord)
        xs = [CR.SVector(x[i], y[i], z[i]) for i in eachindex(x)]
        got = CR.classify_octahedral(xs, 4π .* w)
        @test got !== nothing
        s, θ = got
        @test CR.npoints(s) == length(x)
        # Lebedev's rules are square in the invariants: one unknown per equation
        @test CR.nunknowns(s) == length(CR.invariant_exponents(ord))
        # and the structure reproduces the rule it came from
        xs2, _ = CR.expand(s, θ)
        @test length(xs2) == length(xs)
    end
end

@testset "a family's licence travels automatically" begin
    # Declaring `family_license` is enough: `rule` stamps it into the provenance, so a data
    # package cannot produce rules that look unencumbered by forgetting to set it.
    r = rule(UpstreamLebedev(), Sphere{3}(); degree = 11)
    @test provenance(r).license == CR.family_license(UpstreamLebedev())
    # a family declaring nothing keeps the package's own licence
    @test CR.family_license(LebedevRule()) == ""
    @test provenance(rule(LebedevRule(), Sphere{3}(); degree = 11)).license != ""
    # and the stamp survives the selector's own route
    r2 = rule(Sphere{3}(); degree = 19, copyleft = true)
    @test occursin("GPL-3.0", provenance(r2).license)
end

@testset "large rules are verified on their orbits" begin
    # Degree 29 has 302 points, enough for the grouping to be exercised on every orbit
    # type, and small enough for the general check to finish quickly for comparison.
    r = rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 30)
    fast = check(r)
    full = check(r; use_symmetry = false)
    @test occursin("orbit representatives", fast.basis)
    @test passed(fast) && passed(full)
    @test (fast.exact, fast.sharp, fast.symmetric) == (full.exact, full.sharp, full.symmetric)
end
