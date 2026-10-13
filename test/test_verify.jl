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

@testset "symmetry groups" begin
    # every group, on a rule that has it and on one broken by a moved node
    for (r, g) in ((rule(Simplex{2}(); degree = 6), CR.SymmetricGroup{3}()),
                   (rule(FullySymmetric(), Simplex{3}(); degree = 5), CR.SymmetricGroup{4}()),
                   (rule(Fejer(2), Interval(); degree = 9), CR.Reflection()),
                   (rule(GrundmannMöller(), Simplex{2}(); degree = 7), CR.SymmetricGroup{3}()),
                   (rule(FullySymmetric(), Orthotope{2}(); degree = 9), CR.D4()),
                   (rule(FullySymmetric(), Orthotope{3}(); degree = 7), CR.Oh()),
                   (rule(FullySymmetric(), Orthotope{4}(); degree = 5), CR.B4()),
                   (rule(FullySymmetric(), Disk(); degree = 7), CR.D4()),
                   (rule(FullySymmetric(), Pyramid(); degree = 5), CR.C4v()),
                   (rule(FullySymmetric(), Wedge(); degree = 4), CR.D3h()))
        xs, ws = collect(nodes(r)), collect(weights(r))
        @test provenance(r).symmetry === g
        @test CR.check_invariance(g, xs, ws, 1e-12) === true
        ys = copy(xs)
        ys[end] = ys[end] .+ 1e-6
        @test CR.check_invariance(g, ys, ws, 1e-12) === false
        @test sum(CR.group_orbits(g, xs, ws, 1e-12)[2]) ≈ sum(ws)
    end
    @test CR.check_invariance(CR.NoSymmetry(), [CR.SVector(0.2, 0.3)], [0.5], 1e-12) === nothing
    # a group that does not act on the nodes' dimension: not invariant
    r = rule(Simplex{2}(); degree = 6)
    @test CR.check_invariance(CR.SymmetricGroup{4}(), nodes(r), weights(r), 1e-12) === false
    # orbits whose keys share a first coordinate to rounding are told apart: (½, ½, 0) and
    # (½, ¼, ¼) on the triangle, the first perturbed below and above ½
    pts = [CR.SVector(0.5, 0.5), CR.SVector(0.0, 0.5), CR.SVector(0.5 - 1e-17, 0.0),
           CR.SVector(0.25, 0.25), CR.SVector(0.5, 0.25), CR.SVector(0.25, 0.5)]
    @test CR.check_invariance(CR.SymmetricGroup{3}(), pts, [1.0, 1.0, 1.0, 2.0, 2.0, 2.0], 1e-14) === true
    @test string(CR.SymmetricGroup{3}()) == "S3" && string(CR.Oh()) == "Oh"
end

@testset "in-place BigFloat kernels: outputs are unshared, copies are by value" begin
    setprecision(BigFloat, 128) do
        a, b = BigFloat(1), BigFloat(2)
        v = [a, b, a, b, BigFloat(3)]                 # shared, and not with a neighbour
        CR.unshare!(v, 128)
        @test length(unique(objectid, v)) == 5 && all(x -> precision(x) == 128, v)
        # sharing with numbers outside an output array is the caller's to avoid: copy by value
        dub = CR.DubinerBasis{BigFloat}(3; normalize = true)
        out = zeros(BigFloat, CR.dubiner_length(3))
        CR.dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, big(0.2), big(0.3))
        CR.copy_values!(out, dub.φ)
        kept = copy(out)
        CR.dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, big(0.6), big(0.1))   # writes into dub.φ in place
        @test out == kept && out != dub.φ
    end
end
