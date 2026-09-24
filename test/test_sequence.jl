using CubatureRules, Test
const CR = CubatureRules

@testset "RuleSequence" begin
    s = RuleSequence(ConicalProduct(), Simplex{2}())
    @test length(s) >= 4
    ds = [degree(r) for r in s]
    @test issorted(ds) && allunique(ds)
    @test ds == s.degrees
    f(x) = exp(x[1] * x[2])
    vals = [integrate(f, r) for r in s]
    exact = last(vals)
    @test abs(vals[2] - exact) > abs(vals[end - 1] - exact)         # converging
    @test occursin("RuleSequence", sprint(show, s))
    # an explicit schedule, and BigFloat throughout
    sb = RuleSequence(ConicalProduct(), Simplex{2}(); degrees = [3, 9], digits = 40)
    rs = collect(sb)
    @test [degree(r) for r in rs] == [3, 9]
    @test all(r -> eltype(r) === BigFloat, rs)
    @test_throws ArgumentError RuleSequence(XiaoGimbutas(), Simplex{2}(); degrees = Int[])
end

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

@testset "accuracy-targeted integrate" begin
    f(x) = exp(x[1] * x[2])
    for rtol in (1e-6, 1e-12)
        r = integrate(f, Simplex{2}(); rtol)
        @test r isa IntegrationResult
        @test r.converged
        @test r.error_estimate <= rtol * abs(r.value)
        @test r.neval > 0
    end
    # the estimate is honest: compare with a much finer fixed rule
    ref = integrate(f, rule(Simplex{2}(); degree = 25))
    r = integrate(f, Simplex{2}(); rtol = 1e-8)
    @test abs(r.value - ref) <= 1e-7 * abs(ref)
    # a tolerance it cannot reach reports the last difference and says so
    hard = integrate(x -> exp(x[1]), Simplex{3}(); rtol = 1e-14, maxdegree = 15)
    @test !hard.converged
    @test isfinite(hard.error_estimate)
    # no target at all is an error that explains the alternative
    err = try
        integrate(f, Simplex{2}())
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("rtol", sprint(showerror, err)) && occursin("rule(domain; degree", sprint(showerror, err))
    # an explicit family is honoured
    r = integrate(x -> exp(x), Interval(); rtol = 1e-10, family = Fejer(2))
    @test r.family == "Fejer"
    @test r.converged
end

@testset "a sequence never walks the same rule twice" begin
    # GaussKronrod answers both degree 1 and degree 3 with its 3-point rule. Walking it twice
    # made the successive difference exactly zero, so integrate reported convergence on the
    # first rule it built: exp on [-1,1] came back 6.5e-5 wrong with an error estimate of 0.0.
    seq = RuleSequence(GaussKronrod(), Interval())
    degs = [degree(r) for r in seq]
    @test allunique(degs)
    @test issorted(degs)

    exact = exp(1) - exp(-1)
    res = integrate(x -> exp(x[1]), Interval(); rtol = 1e-10, atol = 1e-14)
    @test res.converged
    @test abs(res.value - exact) <= 1e-10 * exact
    @test res.error_estimate > 0          # a genuine difference, not a rule compared with itself
end
