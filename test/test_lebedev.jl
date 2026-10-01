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
    # Every entry is built and checked in the full sweep (CUBATURERULES_FULL_SWEEP=1); per
    # commit, everything up to degree 17 and the two ends of the 133-201 range, whose rules
    # have 6074 to 13680 points.
    high = filter(e -> e.degree > 17, shipped)
    full = get(ENV, "CUBATURERULES_FULL_SWEEP", "0") == "1"
    ends = isempty(high) ? Int[] : [first(high).degree, last(high).degree]
    for e in shipped
        full || e.degree <= 17 || e.degree in ends || continue
        d = e.degree
        r = rule(LebedevRule(), Sphere{3}(); degree = d)
        v = check(r)
        @test npoints(r) == e.npoints
        @test degree(r) == d
        @test v.exact
        if d <= 17
            @test v.sharp === true
        else
            # in Float64 these rules cannot show that they miss degree d + 1; the miss measured
            # at high precision, recorded in the table, says so instead (see next_error)
            @test v.sharp === nothing && occursin("not resolvable", v.sharp_note)
        end
        @test passed(v)
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

@testset "Lebedev degrees 133 to 201" begin
    high = filter(e -> e.degree > 17, CR.lebedev_entries())
    @test [e.degree for e in high] == collect(133:2:201)
    for e in high
        # square systems: as many unknowns as invariant conditions
        @test CR.nunknowns(e.structure) == length(CR.invariant_exponents(e.degree))
        wmin, dmin = CR.octahedral_margins(e.structure, e.seed)
        @test wmin > 0 && dmin > 0
        # 0.2-1.5% above the parameter-count estimate (d+1)²/3, which is not a proven bound
        @test e.npoints <= 1.02 * (e.degree + 1)^2 / 3
        # certified, and sharp at high precision
        @test e.residual < 1e-14 && e.residual_bits >= 256
        @test e.next_error > 0
    end
    # nothing is tabulated between degrees 17 and 133, so the degree-133 rule answers there
    @test CR.lebedev_entry_for(19).degree == 133
    @test npoints(LebedevRule(), Sphere{3}(), 131) == 6074
    # the selector takes it only where it has fewer points than the product rule
    @test family(rule(Sphere{3}(); degree = 41)) == "SphereProduct"
    @test family(rule(Sphere{3}(); degree = 121)) == "Lebedev"
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

@testset "verification on orbit representatives" begin
    # The test set at degree d is the sorted even exponents of one even degree, which must
    # number exactly the invariants of that degree — 352 at degree 125.
    b = CR.OctahedralTestSet{BigFloat}(125, true)
    @test length.(CR.blocks(b)) == [length(CR.invariant_exponents(125)), 363]

    # The orbit check must agree with the general one — all harmonics at every node, and
    # the all-pairs symmetry check — on every shipped rule up to degree 17, refined to 30
    # digits. Above that the table jumps to degrees 133-201, where refinement takes minutes
    # per rule and the general check evaluates tens of thousands of harmonics at every node.
    for d in 3:2:17
        r = rule(LebedevRule(), Sphere{3}(); degree = d, digits = 30)
        fast, full = check(r), check(r; use_symmetry = false)
        @test occursin("orbit representatives", fast.basis)
        @test (fast.exact, fast.sharp, fast.symmetric, passed(fast)) ==
              (full.exact, full.sharp, full.symmetric, passed(full))
    end

    # And it must fail the rules it should fail.
    r = rule(LebedevRule(), Sphere{3}(); degree = 11)
    relabel(q, c) = QuadratureRule(nodes(q), weights(q), domain(q), c, provenance(q))
    over = check(relabel(r, PolynomialDegree(13)))
    @test !over.exact && !passed(over)
    under = check(relabel(r, PolynomialDegree(9)))
    @test under.exact && under.sharp === false && !passed(under)
    # an even claim is understated by central symmetry alone
    @test check(relabel(r, PolynomialDegree(10))).sharp === false
    # a broken symmetry cannot be grouped into orbits, so the general check takes over
    w = copy(weights(r)); w[1] *= 1.001
    asym = check(QuadratureRule(nodes(r), w, domain(r), PolynomialDegree(11), provenance(r)))
    @test asym.symmetric === false && !passed(asym)
    @test !occursin("orbit representatives", asym.basis)
    @test CR.octahedral_orbits([big.(collect(x)) for x in nodes(r)], big.(w), big(1e-12)) === nothing
end
