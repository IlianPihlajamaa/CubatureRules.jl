using CubatureRules, Test, StaticArrays
const CR = CubatureRules

@testset "Orthotope" begin
    o = Orthotope{2}()
    @test isreference(o)
    @test measure(o) == 4
    @test measure(Orthotope{3}()) == 8
    b = Orthotope((0.0, 1.0), (2.0, 3.0))
    @test !isreference(b)
    @test measure(b) ≈ 4.0
    @test isinterior(SVector(1.0, 2.0), b)
    @test !isinterior(SVector(0.0, 2.0), b)
    @test indomain(SVector(0.0, 2.0), b)
    @test !indomain(SVector(-0.1, 2.0), b)
    @test_throws ArgumentError Orthotope((1.0, 0.0), (0.0, 1.0))
    @test_throws ArgumentError Orthotope((0.0,), (1.0, 2.0))
    @test monomial_moment(o, (2, 4)) == (2 // 3) * (2 // 5)
    @test monomial_moment(o, (1, 2)) == 0
    @test sprint(show, o) == "Orthotope{2}()"
    @test occursin("Orthotope((0.0, 1.0)", sprint(show, b))
end

@testset "tensor-product rules" begin
    for d in (1, 5, 12)
        r = rule(Orthotope{2}(); degree = d)
        @test family(r) == "TensorProduct"
        @test npoints(r) == cld(d + 1, 2)^2
        @test sum(weights(r)) ≈ 4
        v = check(r)
        @test v.exact && v.sharp === true && v.positive && v.interior
    end
    r = rule(Orthotope{3}(); degree = 7, digits = 40)
    @test npoints(r) == 4^3
    @test passed(check(r))
    @test abs(integrate(x -> x[1]^2 * x[2]^4 * x[3]^0, r) - big(2) / 3 * big(2) / 5 * 2) < big(10.0)^-38
    # the claim is the conservative total degree, not the tensor-product space
    @test degree(rule(Orthotope{2}(); degree = 6)) == 7          # 4 points per axis
    @test integrate(x -> x[1]^7 * x[2]^7, rule(Orthotope{2}(); degree = 7)) ≈ 0 atol = 1e-14
    # …though it is in fact exact well beyond the claim, which is why the claim is stated
    @test integrate(x -> x[1]^6 * x[2]^6, rule(Orthotope{2}(); degree = 7)) ≈ (2 / 7)^2 rtol = 1e-14
    # mapped onto a general box
    b = Orthotope((0.0, 1.0), (2.0, 3.0))
    rb = rule(b; degree = 6)
    @test domain(rb) == b
    @test sum(weights(rb)) ≈ measure(b)
    @test passed(check(rb))
    @test integrate(x -> x[1] * x[2], rb) ≈ 2.0 * 4.0 rtol = 1e-13
end

@testset "⊗" begin
    a = rule(Interval(); degree = 5)
    b = rule(Interval(); degree = 9)
    t = a ⊗ b
    @test npoints(t) == npoints(a) * npoints(b)
    @test degree(t) == 5                                  # conservative
    @test domain(t) == Orthotope{2,Float64}(SVector(-1.0, -1.0), SVector(1.0, 1.0))
    @test sum(weights(t)) ≈ 4
    @test passed(check(t))
    @test integrate(x -> x[1]^4 * x[2]^4, t) ≈ (2 / 5)^2 rtol = 1e-14
    # a box rule tensored with an interval rule gives a 3-box
    t3 = t ⊗ rule(Interval(); degree = 3)
    @test CR.dimension(domain(t3)) == 3
    @test degree(t3) == 3
    @test sum(weights(t3)) ≈ 8
    @test_throws ArgumentError rule(Simplex{2}(); degree = 4) ⊗ a
end

@testset "registry on boxes" begin
    a = available(Orthotope{2}(); degree = 9)
    @test !isempty(a)
    @test first(a).family == "TensorProduct(GaussLegendre)"
    @test first(a).npoints == 25
    @test isempty(candidates(TensorProduct, Orthotope((0.0, 0.0), (1.0, 1.0)), PolynomialDegree(3)))
    @test isempty(candidates(TensorProduct, Simplex{2}(), PolynomialDegree(3)))
    # exact rational box rules: a tensor product over Newton–Cotes
    q = rule(Orthotope{2}(); degree = 3, T = Rational{BigInt})
    @test eltype(q) == Rational{BigInt}
    @test sum(weights(q)) == 4
    @test integrate(x -> x[1]^2 * x[2]^2, q) == (2 // 3)^2
    @test check(q).exact
    # explicit construction, including a mixed-family product
    f = TensorProduct((GaussLegendre(), GaussJacobi(1, 0)))
    @test occursin("⊗", CR.describe_family(f))
    @test npoints(f, Orthotope{2}(), 5) == 9
    r = rule(TensorProduct(GaussLegendre(), 2), Orthotope{2}(); degree = 5)
    @test passed(check(r))
    @test occursin("conservative total degree", join(provenance(r).path, " "))
end
