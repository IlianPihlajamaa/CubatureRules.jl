using CubatureRules, Test
const CR = CubatureRules


@testset "embedded rules and error estimates" begin
    e = embedded(GaussKronrod(), Interval(); degree = 11)
    @test e isa EmbeddedRule
    @test npoints(e) == npoints(e.fine)
    @test degree(e) > degree(e.coarse_claim)   # the coarse rule is genuinely coarser
    @test count(!iszero, e.coarse_weights) == CR.kronrod_halves(11)     # the Gauss subset
    exact = exp(1) - exp(-1)
    v = integrate(exp, e)
    @test v ≈ exact rtol = 1e-12
    res = integrate(exp, e; error = true)
    @test res isa IntegrationResult
    @test res.value == v
    @test res.neval == npoints(e)                                       # no extra evaluations
    @test res.error_estimate >= abs(res.value - exact)                  # conservative, as intended
    @test occursin("IntegrationResult", sprint(show, MIME"text/plain"(), res))
    # the coarse rule must really be embedded
    @test_throws ArgumentError EmbeddedRule(rule(Interval(); degree = 5), rule(Interval(); degree = 3))
    # Fejér 2 nests: the 3-point rule's nodes are a subset of the 7-point rule's, so a pair
    # built by hand also estimates its error for free
    fe = EmbeddedRule(rule(Fejer(2), Interval(); degree = 7), rule(Fejer(2), Interval(); degree = 3))
    @test count(!iszero, fe.coarse_weights) == 3
    rf = integrate(exp, fe; error = true)
    @test rf.neval == 7
    @test rf.error_estimate >= abs(rf.value - exact)
    @test rf.value ≈ exact rtol = 1e-5
end

@testset "integrate takes a rule, not a domain" begin
    # There is no tolerance-driven driver; asking for one says what to do instead.
    e = try
        integrate(x -> exp(x), Interval(); rtol = 1e-10)
    catch err
        err
    end
    @test e isa ArgumentError
    msg = sprint(showerror, e)
    @test occursin("rule(", msg)
    @test occursin("QuadGK", msg) && occursin("HCubature", msg)
    @test_throws ArgumentError integrate(x -> x[1], Simplex{2}())
end
