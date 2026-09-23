using CubatureRules, Test
import CubatureRules: ExplicitSeed, MultistartSeed, candidates   # public, not exported
const CR = CubatureRules

# Lebedev's published counts, used here only to document where ours agree and where they do
# not — never as input to the search that produced the shipped seeds.
const LEBEDEV_PUBLISHED = Dict(3 => 6, 5 => 14, 7 => 26, 9 => 38, 11 => 50, 13 => 74, 15 => 86,
                               17 => 110, 19 => 146, 21 => 170, 23 => 194, 25 => 230, 27 => 266,
                               29 => 302)

@testset "Lebedev rules" begin
    shipped = CR.lebedev_entries()
    @test !isempty(shipped)
    @test all(e -> isodd(e.degree), shipped)          # even degrees come free of central symmetry
    @test all(e -> e.status == "ok", shipped)
    for e in shipped
        @test CR.npoints(e.structure) == e.npoints
        @test length(e.seed) == CR.nunknowns(e.structure)
    end
    for e in shipped
        d = e.degree
        r = rule(LebedevRule(), Sphere{3}(); degree = d)
        v = check(r)
        @test npoints(r) == e.npoints
        @test degree(r) == d
        @test v.exact && v.sharp === true
        @test v.positive && v.interior && v.weights_sum_ok
        @test v.symmetric === true                     # the full 48-element group, checked directly
        @test sum(weights(r)) ≈ 4π rtol = 1e-13
        @test all(x -> indomain(x, Sphere{3}()), nodes(r))
        # our count against Lebedev's, where his is known
        pub = get(LEBEDEV_PUBLISHED, d, nothing)
        pub === nothing && continue
        if d == 13
            # his 74-point rule has a negative weight; ours is the smallest positive one
            @test e.npoints == 78
            @test occursin("not positive", e.note)
        else
            @test e.npoints == pub
            @test occursin("matches", e.note)
        end
    end
end

@testset "Lebedev construction and selection" begin
    d = maximum(e -> e.degree, CR.lebedev_entries())
    # the selector prefers it to the product rule wherever it exists, on node count
    r = rule(Sphere{3}(); degree = 11)
    @test family(r) == "Lebedev"
    @test npoints(r) < CR.npoints(SphereProduct(), Sphere{3}(), 11)
    @test occursin("selected by rule()", provenance(r).selection)
    # an even-degree request is answered by the odd rule above it, for free
    re = rule(LebedevRule(), Sphere{3}(); degree = 10)
    @test degree(re) == 11
    # arbitrary precision
    for digits in (30, 60)
        rb = rule(LebedevRule(), Sphere{3}(); degree = 9, digits)
        @test passed(check(rb))
        dev = CR.with_bits(() -> abs(sum(weights(rb)) - 4 * BigFloat(π)), CR.digits_to_bits(digits) + 64)
        @test dev < big(10.0)^(-digits + 2)
    end
    @test_throws NoRuleError rule(LebedevRule(), Sphere{3}(); degree = 5, T = Rational{BigInt})
    @test_throws NoRuleError rule(LebedevRule(), Sphere{3}(); degree = d + 2)
    @test isempty(candidates(LebedevRule, Sphere{2}(), PolynomialDegree(5)))
    @test isempty(candidates(LebedevRule, Simplex{2}(), PolynomialDegree(5)))
    # the provenance records the group and the licence of the seeds
    r5 = rule(LebedevRule(), Sphere{3}(); degree = 5)
    @test provenance(r5).symmetry === :Oh
    @test occursin("generated in-house", provenance(r5).license)
    @test occursin("invariants", certificate(r5).equations)
    @test certificate(r5).residual < 1e-14      # a Float64 rule: the residual is at its roundoff
    @test occursin("Lebedev1976", cite(r5))
end

@testset "Lebedev seeds are recoverable" begin
    # the shipped numbers are a convenience, not an input: the rule can be found again from
    # its orbit structure alone
    r = rule(LebedevRule(), Sphere{3}(); degree = 9, seed = MultistartSeed(; nstarts = 256))
    v = check(r)
    @test v.exact && v.positive && v.symmetric === true
    @test npoints(r) == npoints(rule(LebedevRule(), Sphere{3}(); degree = 9))
    @test occursin("multistart", lowercase(provenance(r).seed_source))
    # an explicit seed is accepted and recorded as such
    e = CR.lebedev_entry_for(9)
    r2 = rule(LebedevRule(), Sphere{3}(); degree = 9, seed = ExplicitSeed(e.seed; source = "the shipped table"))
    @test passed(check(r2))
    @test occursin("the shipped table", provenance(r2).seed_source)
    @test_throws ArgumentError rule(LebedevRule(), Sphere{3}(); degree = 9, seed = ExplicitSeed([0.1]))
end

@testset "octahedral symmetry check" begin
    r = rule(LebedevRule(), Sphere{3}(); degree = 5)
    xs, ws = collect(nodes(r)), collect(weights(r))
    @test CR.check_octahedral_symmetry(xs, ws, 1e-12)
    # break one weight: the group orbit no longer carries a single weight
    ws2 = copy(ws); ws2[1] *= 1.5
    @test !CR.check_octahedral_symmetry(xs, ws2, 1e-12)
    # move one node off its orbit
    xs2 = copy(xs); xs2[1] = CR.SVector(0.0, 0.6, 0.8)
    @test !CR.check_octahedral_symmetry(xs2, ws, 1e-12)
    # a rule on the circle is not in the group's dimension
    @test !CR.check_octahedral_symmetry([CR.SVector(1.0, 0.0)], [1.0], 1e-12)
end
