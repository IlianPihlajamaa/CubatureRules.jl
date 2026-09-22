using CubatureRules, Test, StaticArrays, LinearAlgebra
const CR = CubatureRules

@testset "claims" begin
    @test degree(PolynomialDegree(7)) == 7
    @test_throws ArgumentError PolynomialDegree(-1)
    @test_throws ArgumentError degree(NoClaim())
    @test_throws ArgumentError degree(SpanOf([sin, cos]))
    @test SpanOf([1, 2]) == SpanOf([1, 2])
end

@testset "domains" begin
    s = Simplex{2}()
    @test isreference(s)
    @test measure(s) == 1 // 2
    @test measure(Simplex{3}()) == 1 // 6
    t = Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 3.0))
    @test !isreference(t)
    @test measure(t) ≈ 3.0
    @test barycentric(t, SVector(1.0, 1.0)) ≈ SVector(1 - 0.5 - 1 / 3, 0.5, 1 / 3)
    @test cartesian(t, SVector(0.2, 0.3, 0.5)) ≈ SVector(0.6, 1.5)
    @test isinterior(SVector(0.2, 0.2), s)
    @test !isinterior(SVector(0.0, 0.2), s)
    @test indomain(SVector(0.0, 0.2), s)
    @test !indomain(SVector(0.6, 0.6), s)
    @test indomain(SVector(0.6, 0.41), s; tol = 0.02)
    @test measure(Interval(0, 3)) == 3
    @test isreference(Interval())
    @test_throws ArgumentError Interval(1, 0)
    @test_throws CR.NotYetImplemented Polytope{3}()
    @test_throws CR.NotYetImplemented Wedge()
    w = WeightedDomain(Interval(), JacobiWeight(1, 0))
    @test measure(w) == 2
end

@testset "exact moments (Dirichlet)" begin
    @test monomial_moment(Simplex{2}(), (0, 0)) == 1 // 2
    @test monomial_moment(Simplex{2}(), (1, 0)) == 1 // 6
    @test monomial_moment(Simplex{2}(), (2, 3)) == factorial(2) * factorial(3) // factorial(big(7))
    @test monomial_moment(Simplex{3}(), (1, 1, 1)) == 1 // factorial(big(6))
    @test monomial_moment(Simplex{4}(), (0, 0, 0, 0)) == 1 // 24
    @test barycentric_moment((1, 1, 1)) == 1 // 120
    @test monomial_moment(Interval(), (4,)) == 2 // 5
    @test monomial_moment(Interval(), (3,)) == 0
end

@testset "orthonormal Dubiner basis" begin
    n = 8
    r = rule(ConicalProduct(), Simplex{2}(); degree = 2n + 1, digits = 40)
    L = CR.dubiner_length(n)
    ws = CR.DubinerWorkspace{BigFloat}(n)
    φ = zeros(BigFloat, L)
    G = zeros(BigFloat, L, L)
    for (x, w) in zip(r.nodes, r.weights)
        CR.dubiner!(φ, ws, x[1], x[2])
        G .+= w .* φ .* φ'
    end
    @test maximum(abs, G - I) < 1e-35
    # gradients against finite differences
    ws64 = CR.DubinerWorkspace{Float64}(n)
    v, gx, gy = zeros(L), zeros(L), zeros(L)
    CR.dubiner!(v, gx, gy, ws64, 0.3, 0.25)
    h = 1e-6
    vp, vm = zeros(L), zeros(L)
    CR.dubiner!(vp, ws64, 0.3 + h, 0.25); CR.dubiner!(vm, ws64, 0.3 - h, 0.25)
    @test maximum(abs, (vp - vm) / 2h - gx) < 1e-6
    CR.dubiner!(vp, ws64, 0.3, 0.25 + h); CR.dubiner!(vm, ws64, 0.3, 0.25 - h)
    @test maximum(abs, (vp - vm) / 2h - gy) < 1e-6
end

@testset "rule records, static form, equality" begin
    r = rule(Simplex{2}(); degree = 4)
    @test npoints(r) == length(nodes(r)) == length(weights(r))
    @test degree(r) >= 4
    @test domain(r) == Simplex{2}()
    @test exactness(r) isa PolynomialDegree
    @test provenance(r) isa Provenance
    @test certificate(r) isa Certificate
    s = static(r)
    @test s isa StaticQuadratureRule
    @test isbits(s)
    @test npoints(s) == npoints(r)
    @test family(s) == family(r)
    # 1D specialisation stores scalars
    l = rule(Interval(); degree = 7)
    @test nodes(l) isa Vector{Float64}
    @test isbits(static(l))
    # content equality and hashing: same construction twice is the same rule
    r2 = rule(Simplex{2}(); degree = 4)
    @test r == r2
    @test hash(r) == hash(r2)
    @test isapprox(r, r2)
    @test Dict(r => 1)[r2] == 1
    @test rule_hash(r) == rule_hash(r2)
    @test r != rule(Simplex{2}(); degree = 5)
    @test_throws DimensionMismatch QuadratureRule([SVector(0.1, 0.1)], [0.1, 0.2], Simplex{2}(), NoClaim(),
                                                  Provenance(family = "x", derivation = Derived()))
end

@testset "precision resolution" begin
    @test CR.resolve_precision(nothing, nothing) == (Float64, 53)
    @test CR.resolve_precision(nothing, 100)[1] === BigFloat
    @test CR.resolve_precision(nothing, 100)[2] >= 333
    @test CR.resolve_precision(Float32, nothing) == (Float32, 24)
    @test setprecision(() -> CR.resolve_precision(BigFloat, nothing), BigFloat, 300) == (BigFloat, 300)
    @test_throws ArgumentError CR.resolve_precision(Float64, 30)
    # the ambient precision is read once, at the call; generators never depend on it
    r = setprecision(() -> rule(Interval(); degree = 5, T = BigFloat), BigFloat, 150)
    @test precision(first(weights(r))) == 150
    r2 = rule(Interval(); degree = 5, digits = 30)
    @test precision(first(weights(r2))) == CR.digits_to_bits(30)
end

@testset "orthonormal tetrahedral Dubiner basis" begin
    n = 6
    r = rule(ConicalProduct(), Simplex{3}(); degree = 2n + 1, digits = 40)
    b = CR.SimplexBasis{3,BigFloat}(n)
    L = CR.basis_length(b)
    @test L == CR.tet_length(n) == binomial(n + 3, 3)
    G = zeros(BigFloat, L, L)
    for (x, w) in zip(r.nodes, r.weights)
        φ, _ = CR.evaluate!(b, x; gradient = false)
        G .+= w .* φ .* φ'
    end
    @test maximum(abs, G - I) < 1e-35
    # gradients against finite differences
    b64 = CR.SimplexBasis{3,Float64}(n)
    x0 = [0.2, 0.15, 0.3]
    _, Gr = CR.evaluate!(b64, x0)
    Gr = copy(Gr)
    h = 1e-6
    for j in 1:3
        xp = copy(x0); xp[j] += h
        xm = copy(x0); xm[j] -= h
        fp = copy(CR.evaluate!(b64, xp; gradient = false)[1])
        fm = copy(CR.evaluate!(b64, xm; gradient = false)[1])
        @test maximum(abs, (fp - fm) / 2h - Gr[:, j]) < 1e-6
    end
    # degree blocks are contiguous and complete
    @test [length(CR.degree_block(3, k)) for k in 0:4] == [(k + 1) * (k + 2) ÷ 2 for k in 0:4]
    # exact (unnormalised) evaluation stays rational
    br = CR.SimplexBasis{3,Rational{BigInt}}(4; normalize = false)
    φ, _ = CR.evaluate!(br, Rational{BigInt}[1 // 5, 1 // 7, 1 // 3])
    @test eltype(φ) == Rational{BigInt}
end
