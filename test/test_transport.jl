# Claim preservation (PLAN §2.3): every transport states what happens to the claim.
using CubatureRules, Test, StaticArrays, LinearAlgebra
const CR = CubatureRules

@testset "map_to (affine) preserves the claim" begin
    r = rule(Simplex{2}(); degree = 8)
    t = Simplex((1.0, 0.0), (3.0, 1.0), (0.5, 2.0))
    m = map_to(r, t)
    @test exactness(m) == exactness(r)
    @test domain(m) == t
    @test sum(weights(m)) ≈ measure(t)
    @test passed(check(m))
    @test occursin("map_to", last(provenance(m).path))
    # exact rational rules stay exact under rational maps
    g = rule(GrundmannMöller(), Simplex{2}(); degree = 5, T = Rational{BigInt})
    gm = map_to(g, Simplex((0, 0), (2, 0), (0, 3)))
    @test eltype(gm) == Rational{BigInt}
    @test check(gm).exact
    li = map_to(rule(Interval(); degree = 7), Interval(0.0, 5.0))
    @test exactness(li) == PolynomialDegree(7)
    @test passed(check(li))
    @test_throws ArgumentError map_to(r, Interval())
end

@testset "subdivide preserves the claim" begin
    r = rule(Simplex{2}(); degree = 4)
    s = subdivide(r, 3)
    @test npoints(s) == 9npoints(r)
    @test exactness(s) == PolynomialDegree(4)
    @test sum(weights(s)) ≈ 0.5
    @test check(s).exact
    t = Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
    st = subdivide(r, t, 2)
    @test domain(st) == t
    @test check(st).exact
    li = subdivide(rule(Interval(); degree = 5), Interval(0.0, 1.0), 4)
    @test npoints(li) == 12
    @test check(li).exact
    # composite rules converge on non-polynomial integrands
    f(x) = sqrt(x[1] + x[2])
    e1 = abs(integrate(f, subdivide(r, 2)) - integrate(f, subdivide(r, 32)))
    e2 = abs(integrate(f, subdivide(r, 8)) - integrate(f, subdivide(r, 32)))
    @test e2 < e1
end

@testset "transform and duffy destroy the claim" begin
    r = rule(Simplex{2}(); degree = 6)
    φ(x) = SVector(x[1]^2, x[2])
    Jφ(x) = [2x[1] 0; 0 1]
    tr = transform(r, φ, Jφ)
    @test exactness(tr) isa NoClaim
    @test_throws ArgumentError degree(tr)
    @test_throws ArgumentError verify(tr)
    @test occursin("NoClaim", last(provenance(tr).path))
    # a caller may assert a claim; it is recorded as an assertion
    ta = transform(r, x -> x, x -> 1.0; claim = PolynomialDegree(6))
    @test exactness(ta) == PolynomialDegree(6)
    @test occursin("ASSERTED", last(provenance(ta).path))
    d = duffy(r)
    @test exactness(d) isa NoClaim
    @test provenance(d).symmetry === :none
    # radial grading removes the 1/|x| blow-up at the vertex (the angular dependence stays)
    sing(x) = 1 / sqrt(x[1]^2 + x[2]^2)
    exact_val = sqrt(2) * log(1 + sqrt(2))         # ∫_T 1/|x| over the reference triangle
    # a convergence sweep is the only check available, and it is empirical
    seq = [duffy(rule(ConicalProduct(), Simplex{2}(); degree = d)) for d in (3, 7, 15, 31)]
    # grading fixes the radial blow-up but not the angular dependence, so what is left
    # converges algebraically: 2.5e-2, 3.7e-3, 3.3e-4, 2.5e-5 over these four degrees. The
    # target has to say that, rather than the √eps a smooth integrand would reach.
    v = CR.verify_convergence(seq, sing, exact_val; rtol = 1e-4)
    @test v.empirical && v.method === :convergence_sweep
    @test v.exact
end

@testset "subdivide on tetrahedra" begin
    r = rule(Simplex{3}(); degree = 4)
    s = subdivide(r, 2)
    @test npoints(s) == 8npoints(r)
    @test exactness(s) == exactness(r)
    @test sum(weights(s)) ≈ 1 / 6
    @test check(s).exact
    @test_throws ArgumentError subdivide(r, 3)
    # bisection keeps the cells non-degenerate and space-filling
    cells = CR.bisect_tetrahedron(Simplex{3,Float64,4}(CR.SVector{4}(map(v -> CR.SVector{3,Float64}(v),
                                                                        Simplex{3}().vertices))))
    @test length(cells) == 8
    @test sum(measure, cells) ≈ 1 / 6
    @test all(c -> measure(c) ≈ 1 / 48, cells)
    t = Simplex((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 3.0))
    st = subdivide(r, t, 2)
    @test sum(weights(st)) ≈ measure(t)
    @test check(st).exact
end
