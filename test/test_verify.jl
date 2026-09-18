using CubatureRules, Test, StaticArrays
const CR = CubatureRules

@testset "a correct rule passes" begin
    v = check(rule(Simplex{2}(); degree = 11))
    @test passed(v)
    @test v.method === :exact_integration && !v.empirical
    @test v.precision_bits == 128
end

@testset "verification catches what certification cannot" begin
    r = rule(Simplex{2}(); degree = 10)
    # a transcription error in one weight: exactness fails
    w = copy(weights(r)); w[3] *= 1 + 1e-9
    bad = QuadratureRule(nodes(r), w, domain(r), exactness(r), provenance(r), certificate(r))
    vb = check(bad)
    @test !vb.exact
    @test !passed(vb)
    @test vb.symmetric === false          # the perturbed weight also breaks the orbit
    # an overstated claim fails exactness
    over = QuadratureRule(nodes(r), weights(r), domain(r), PolynomialDegree(11), provenance(r), certificate(r))
    @test !check(over).exact
    # an understated claim fails sharpness
    under = QuadratureRule(nodes(r), weights(r), domain(r), PolynomialDegree(9), provenance(r), certificate(r))
    vu = check(under)
    @test vu.exact
    @test vu.sharp === false
    @test !passed(vu)
    # testing at a user-chosen degree skips sharpness
    @test check(r; degree = 7).sharp === nothing
    @test check(r; degree = 7).exact
end

@testset "positivity and interiority are reported" begin
    g = rule(GrundmannMöller(), Simplex{2}(); degree = 7)
    v = check(g)
    @test !v.positive
    @test v.interior
    @test passed(v)                      # negative weights are not a claim violation
end

@testset "higher-precision checks" begin
    r = rule(Simplex{2}(); degree = 14, digits = 60)
    v = check(r)
    @test v.precision_bits == 2precision(first(weights(r)))
    @test passed(v)
    @test check(r; bits = 1000).precision_bits == 1000
end

@testset "SpanOf claims" begin
    # the 2-point Gauss rule is exact on {1, x, x², x³}; claim a span of two of them
    l = rule(Interval(); degree = 3)
    s = QuadratureRule(nodes(l), weights(l), domain(l), SpanOf([x -> one(x), x -> x^2]), provenance(l))
    @test verify(s; integrals = [2, 2 // 3]).exact
    @test !verify(s; integrals = [2, 1]).exact
    @test_throws ArgumentError verify(s)
end

@testset "@test_exact" begin
    r = rule(Simplex{2}(); degree = 9)
    @test_exact r 9
    @test_exact r 5
end
