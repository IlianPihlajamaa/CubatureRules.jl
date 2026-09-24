using CubatureRules, Test
const CR = CubatureRules

@testset "Gauss–Patterson" begin
    # the 1- and 3-point levels are the midpoint rule and 3-point Gauss
    @test nodes(rule(GaussPatterson(), Interval(); degree = 1)) == [0.0]
    r3 = rule(GaussPatterson(), Interval(); degree = 5)
    @test nodes(r3) ≈ [-sqrt(0.6), 0, sqrt(0.6)]
    @test weights(r3) ≈ [5 / 9, 8 / 9, 5 / 9]

    # the published 7-point rule (Patterson 1968)
    r7 = rule(GaussPatterson(), Interval(); npoints = 7)
    @test degree(r7) == 11
    @test nodes(r7)[5:7] ≈ [0.4342437493468025580020715, 0.7745966692414833770358531, 0.9604912687080202834235071] atol = 1e-15
    @test weights(r7)[4:7] ≈ [0.4509165386584741423451, 0.4013974147759622229051,
                              0.2684880898683334407286, 0.1046562260264672651938] rtol = 1e-14

    # levels: 2^(k+1) - 1 points, degree 3·2^k - 1
    @test [npoints(rule(GaussPatterson(), Interval(); degree = d)) for d in (0, 1, 2, 5, 6, 11, 12, 23, 24, 47)] ==
          [1, 1, 3, 3, 7, 7, 15, 15, 31, 31]
    @test_throws ArgumentError rule(GaussPatterson(), Interval(); npoints = 8)
    @test_throws NoRuleError rule(GaussPatterson(), Interval(); degree = 768)

    # each level contains the one before, has positive weights and verifies
    prev = nothing
    for n in (3, 7, 15, 31, 63)
        r = rule(GaussPatterson(), Interval(); npoints = n)
        @test CR.passed(verify(r))
        @test all(>(0), weights(r))
        prev === nothing || @test all(p -> any(q -> abs(p - q) <= 4eps(), nodes(r)), nodes(prev))
        prev = r
    end
    @test CR.passed(verify(rule(GaussPatterson(), Interval(); npoints = 63, digits = 50)))

    # From 127 points on a rule misses the next degree by less than Float64 rounding. That is
    # reported as unresolvable at this precision, not as an understated claim, and at enough
    # digits the miss shows.
    r127 = rule(GaussPatterson(), Interval(); npoints = 127)
    v = verify(r127)
    @test CR.passed(v)
    @test v.sharp === nothing && occursin("not resolvable", v.sharp_note)
    @test 1e-22 < certificate(r127).next_error < 1e-18
    @test occursin("not resolvable", sprint(show, MIME"text/plain"(), v))
    @test verify(rule(GaussPatterson(), Interval(); npoints = 127, digits = 40)).sharp === true
    # a rule without the record still fails when it is exact one degree too high
    g = rule(GaussLegendre(), Interval(); npoints = 5)
    understated = QuadratureRule(nodes(g), weights(g), Interval(), PolynomialDegree(8), provenance(g), certificate(g))
    @test verify(understated).sharp === false

    # the embedded pair: the level asked for and the one below, from one set of values
    e = embedded(GaussPatterson(), Interval(); degree = 20)
    @test npoints(e) == 15 && count(!iszero, e.coarse_weights) == 7
    res = integrate(x -> 1 / (1 + 25x^2), e; error = true)
    @test res.error_estimate >= abs(res.value - 2atan(5) / 5)
    @test_throws ArgumentError embedded(GaussPatterson(), Interval(); degree = 1)

    # mapped intervals, and the selector still prefers Gauss for a plain degree request
    @test integrate(x -> x^11, rule(GaussPatterson(), Interval(0, 2); degree = 11)) ≈ 2^12 / 12
    @test family(rule(Interval(); degree = 11)) != "GaussPatterson"
end
