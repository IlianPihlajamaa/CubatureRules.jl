# Interop with Lebedev.jl, which is GPL-3. Loading it here is *use*, not distribution: no
# part of it is redistributed by this package, and the test environment is not shipped.
#
# Imported under an alias on purpose. This package exports a type named `Lebedev`, so
# `using Lebedev` would collide with it; the extension loads either way.
using CubatureRules, Test
import Lebedev as LebedevJL
const CR = CubatureRules

@testset "one family, two sources of seeds" begin
    @test UpstreamLebedev === CR.Lebedev{CR.LebedevJLSeeds}
    @test Lebedev() isa Lebedev            # the in-house default
    @test Lebedev() === CR.Lebedev{CR.InHouseSeeds}()
    @test CR.family_name(UpstreamLebedev()) == CR.family_name(Lebedev()) == "Lebedev"
    # they are told apart by the description, which carries the licence
    @test CR.describe_family(Lebedev()) == "Lebedev"
    @test occursin("GPL-3", CR.describe_family(UpstreamLebedev()))
    @test CR.family_license(Lebedev()) == ""
    @test occursin("GPL-3.0", CR.family_license(UpstreamLebedev()))
    @test CR.missing_dependency(UpstreamLebedev()) === nothing      # the extension is loaded
end

@testset "available lists both; rule takes only one" begin
    names = [row.family for row in available(Sphere{3}(); degree = 11)]
    @test "Lebedev" in names                       # ours
    @test "Lebedev (Lebedev.jl, GPL-3)" in names   # theirs, with the terms in the name
    # the selector will not take theirs on its own
    @test CR.selectable(Lebedev()) === true
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
        x, y, z, w = LebedevJL.lebedev_by_order(ord)
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
