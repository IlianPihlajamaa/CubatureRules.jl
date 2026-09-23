using CubatureRules, Test
import CubatureRules: ExponentialWeight, candidates, isreference   # public, not exported
const CR = CubatureRules

@testset "unbounded domains" begin
    @test isreference(HalfLine()) && isreference(RealLine())
    @test measure(HalfLine()) == Inf && measure(RealLine()) == Inf
    @test indomain(0.0, HalfLine()) && !indomain(-1.0, HalfLine())
    @test isinterior(1.0, HalfLine()) && !isinterior(0.0, HalfLine())
    @test isinterior(-1e9, RealLine())
    @test sprint(show, HalfLine()) == "HalfLine()"
    # weighted, the measure is the weight's mass — finite
    @test measure(LaguerreRay()) == 1
    @test measure(LaguerreRay(2)) == 2                       # Γ(3)
    @test measure(HermiteLine()) ≈ sqrt(π)
    @test LaguerreRay().weight(2.0) ≈ exp(-2)
    @test HermiteLine().weight(2.0) ≈ exp(-4)
    @test_throws ArgumentError ExponentialWeight(-2)
    @test_throws ArgumentError GaussLaguerre(-1)
end

@testset "Gauss–Laguerre and Gauss–Hermite" begin
    for (dom, fam) in ((LaguerreRay(), GaussLaguerre()), (LaguerreRay(2), GaussLaguerre(2)),
                       (HermiteLine(), GaussHermite()))
        for d in (1, 3, 9, 21)
            r = rule(dom; degree = d)
            v = check(r)
            @test family(r) == CR.family_name(fam)
            @test domain(r) == dom
            @test npoints(r) == cld(d + 1, 2)
            @test degree(r) == 2npoints(r) - 1
            @test sum(weights(r)) ≈ Float64(measure(dom)) rtol = 1e-13
            @test v.exact && v.sharp === true && v.positive && v.interior
        end
    end
    # the classical moments
    rl = rule(LaguerreRay(); degree = 11)
    for k in 0:5
        @test integrate(x -> x^k, rl) ≈ factorial(k) rtol = 1e-12
    end
    rh = rule(HermiteLine(); degree = 11)
    @test integrate(x -> x^2, rh) ≈ sqrt(π) / 2 rtol = 1e-12
    @test integrate(x -> x^3, rh) ≈ 0 atol = 1e-12          # odd moments vanish
    @test check(rh).symmetric === true                       # Hermite nodes are symmetric
    @test check(rl).symmetric === nothing                    # Laguerre claims none
    # arbitrary precision
    for digits in (30, 60)
        r = rule(HermiteLine(); degree = 15, digits)
        @test passed(check(r))
        @test abs(sum(weights(r)) - sqrt(big(π))) < big(10.0)^(-digits + 2)
        rlag = rule(LaguerreRay(); degree = 15, digits)
        @test passed(check(rlag))
        @test abs(integrate(x -> x^7, rlag) - factorial(big(7))) < big(10.0)^(-digits + 6)
    end
    # requesting by point count, and the selector on these domains
    r = rule(GaussHermite(), HermiteLine(); npoints = 12)
    @test npoints(r) == 12 && degree(r) == 23
    a = available(LaguerreRay(); degree = 9)
    @test first(a).family == "GaussLaguerre"
    @test isempty(candidates(GaussHermite, Interval(), PolynomialDegree(3)))
    @test isempty(candidates(GaussLaguerre, HermiteLine(), PolynomialDegree(3)))
    @test_throws NoRuleError rule(HermiteLine(); degree = 3, T = Rational{BigInt})
end

@testset "exp-sinh and sinh-sinh" begin
    for (fam, dom, lvl) in ((ExpSinh, HalfLine(), 2:5), (SinhSinh, RealLine(), 2:5))
        for m in lvl
            r = rule(fam(m), dom)
            @test exactness(r) isa NoClaim
            @test_throws ArgumentError degree(r)
            @test_throws ArgumentError verify(r)              # nothing to verify by exactness
            @test all(>(0), weights(r))
            @test all(x -> isinterior(x, dom), nodes(r))
            @test issorted(nodes(r))
            @test npoints(r) == npoints(fam(m), dom)
        end
        # never offered for a degree request, and asked for by name instead
        @test isempty(candidates(fam, dom, PolynomialDegree(5)))
        @test isempty(available(dom; degree = 5))
        @test_throws ArgumentError fam(0)
        @test_throws NoRuleError rule(fam(3), dom; T = Rational{BigInt})
    end
    # the model integrands: decay at infinity, integrable singularity at the origin
    es = rule(ExpSinh(5), HalfLine())
    @test family(es) == "ExpSinh"
    @test integrate(x -> exp(-x), es) ≈ 1 atol = 1e-15
    @test integrate(x -> exp(-x) / sqrt(x), es) ≈ sqrt(π) atol = 1e-14
    @test integrate(x -> 1 / (1 + x^2), es) ≈ π / 2 atol = 1e-14
    ss = rule(SinhSinh(5), RealLine())
    @test integrate(x -> exp(-x^2), ss) ≈ sqrt(π) atol = 1e-14
    @test integrate(x -> 1 / (1 + x^2), ss) ≈ π atol = 1e-14
    # the truncation range is set by the output precision, and widens with it
    @test maximum(nodes(es)) > 1e30 && minimum(nodes(es)) < 1e-30
    @test maximum(nodes(rule(ExpSinh(5), HalfLine(); digits = 40))) > big(10.0)^70
    # convergence sweeps, flagged empirical
    seq = [rule(ExpSinh(m), HalfLine()) for m in 2:6]
    v = CR.verify_convergence(seq, x -> exp(-x) / sqrt(x), sqrt(π))
    @test passed(v) && v.empirical && v.method === :convergence_sweep
    @test !passed(CR.verify_convergence(seq, x -> exp(-x) / sqrt(x), 2.0))
    @test passed(CR.verify_convergence([rule(SinhSinh(m), RealLine()) for m in 2:6],
                                       x -> exp(-x^2), sqrt(π)))
    # consecutive levels nest, over sixty orders of magnitude of node
    e = EmbeddedRule(rule(ExpSinh(5), HalfLine()), rule(ExpSinh(4), HalfLine()))
    res = integrate(x -> exp(-x) / sqrt(x), e; error = true)
    @test res.neval == npoints(rule(ExpSinh(5), HalfLine()))
    @test res.degree == -1
    @test res.value ≈ sqrt(π) atol = 1e-14
    @test res.error_estimate < 1e-10
    # arbitrary precision
    rb = rule(SinhSinh(7), RealLine(); digits = 40)
    @test abs(integrate(x -> exp(-x^2), rb) - sqrt(big(π))) < big(10.0)^-38
end

@testset "tolerance-driven integration on unbounded domains" begin
    # nothing here answers a degree request, so the sweep walks levels instead
    r = integrate(x -> exp(-x), HalfLine(); rtol = 1e-10)
    @test r.converged && r.degree == -1 && r.family == "ExpSinh"
    @test r.value ≈ 1 atol = 1e-14
    r = integrate(x -> exp(-x^2), RealLine(); rtol = 1e-10)
    @test r.converged && r.family == "SinhSinh"
    @test r.value ≈ sqrt(π) atol = 1e-13
    r = integrate(x -> exp(-x) / sqrt(x), HalfLine(); rtol = 1e-12)
    @test r.converged
    @test r.value ≈ sqrt(π) atol = 1e-13
    # a weighted domain still has a degree-based family, and keeps using it
    r = integrate(x -> x^5, LaguerreRay(); rtol = 1e-12)
    @test r.converged && r.family == "GaussLaguerre"
    @test r.value ≈ 120 rtol = 1e-13
    # a level family can also be named explicitly, including on an interval
    r = integrate(x -> 1 / sqrt(1 - x^2), Interval(); rtol = 1e-8, family = TanhSinh())
    @test r.converged
    @test r.value ≈ π atol = 1e-7
    s = CR.LevelSequence(ExpSinh(), HalfLine(); levels = 2:5)
    @test length(s) == 4 && eltype(collect(s)) <: QuadratureRule
    @test occursin("levels [2, 3, 4, 5]", sprint(show, s))
    @test_throws ArgumentError CR.LevelSequence(GaussLegendre(), Interval())
    # the selector explains itself when asked for a degree it cannot give
    msg = try rule(HalfLine(); degree = 5) catch e; sprint(showerror, e) end
    @test occursin("ExpSinh", msg) && occursin("rule(ExpSinh(4), HalfLine())", msg)
    msg = try rule(RealLine()) catch e; sprint(showerror, e) end
    @test occursin("SinhSinh", msg)
    @test CR.home_domain(TanhSinh()) == Interval()
    @test CR.home_domain(GaussLegendre()) === nothing
end
