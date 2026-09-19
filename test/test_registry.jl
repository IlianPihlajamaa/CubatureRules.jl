using CubatureRules, Test
const CR = CubatureRules

@testset "discovery" begin
    fams = families()
    @test XiaoGimbutas in fams
    @test GrundmannMöller in fams
    @test any(F -> F <: GaussJacobi, fams)
    @test any(F -> F <: ConicalProduct, fams)
    @test issorted(fams; by = string)             # explicit, deterministic order
    @test !any(F -> F <: CombinatorFamily, CR.leaf_families())
end

@testset "candidates" begin
    @test isempty(candidates(XiaoGimbutas, Interval(), PolynomialDegree(3)))
    @test isempty(candidates(XiaoGimbutas, Simplex{3}(), PolynomialDegree(3)))
    @test isempty(candidates(XiaoGimbutas, Simplex{2}(), PolynomialDegree(60)))
    @test candidates(XiaoGimbutas, Simplex{2}(), PolynomialDegree(20)) == [XiaoGimbutas()]
    # combinators recurse into the 1D leaves
    cp = candidates(ConicalProduct, Simplex{2}(), PolynomialDegree(5))
    @test !isempty(cp) && all(c -> c.inner isa GaussJacobi, cp)
    @test npoints(XiaoGimbutas(), Simplex{2}(), 17) == 60
    @test npoints(ConicalProduct(), Simplex{3}(), 5) == 27
    @test npoints(GrundmannMöller(), Simplex{2}(), 3) == 4
    @test properties(GrundmannMöller(), Simplex{2}(), 5).positive == false
    @test derivation(XiaoGimbutas()) == Seeded()
    @test derivation(GrundmannMöller()) == Derived()
    @test cost_estimate(XiaoGimbutas(), Simplex{2}(), 10, Float64) > 0
end

@testset "available ranks by node count (exit criterion)" begin
    a = available(Simplex{2}(); degree = 17)
    names = [row.family for row in a]
    @test names[1] == "XiaoGimbutas"
    @test "ConicalProduct(GaussJacobi)" in names
    @test "GrundmannMöller" in names
    @test issorted([row.npoints for row in a])
    @test findfirst(==("ConicalProduct(GaussJacobi)"), names) < findfirst(==("GrundmannMöller"), names)
    # filters
    ap = available(Simplex{2}(); degree = 17, positive = true)
    @test !("GrundmannMöller" in [row.family for row in ap])
    # without degree: families and ranges
    a0 = available(Simplex{2}())
    @test length(a0) >= 3
    @test any(row -> row.family == "XiaoGimbutas" && first(row.degrees) == 0 && last(row.degrees) >= 20, a0)
    c = compare(Simplex{2}(), 17)
    @test length(c) >= 3
    @test any(row -> row.exact_rational, c)
end

@testset "rule selects the best candidate and records why" begin
    r = rule(Simplex{2}(); degree = 17)
    @test family(r) == "XiaoGimbutas"
    @test occursin("XiaoGimbutas (60 points) < ConicalProduct", provenance(r).selection)
    @test family(rule(Simplex{2}(); degree = 60)) == "ConicalProduct"      # beyond the seeded range
    @test family(rule(Simplex{2}(); degree = 5, T = Rational{BigInt})) == "GrundmannMöller"
    @test family(rule(Simplex{3}(); degree = 5)) == "FullySymmetric"          # 14 points beats 15 and 27
    @test family(rule(Simplex{2}(); degree = 5, family = ConicalProduct())) == "ConicalProduct"
    @test occursin("explicitly", provenance(rule(ConicalProduct(), Simplex{2}(); degree = 5)).selection)
    # a mapped domain gets the reference rule, mapped affinely
    t = Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
    rt = rule(t; degree = 6)
    @test domain(rt) == t
    @test sum(weights(rt)) ≈ 2.0
end

@testset "diagnostics" begin
    e = try
        rule(Simplex{2}())
    catch err
        err
    end
    @test e isa NoRuleError
    msg = sprint(showerror, e)
    @test occursin("needs a `degree`", msg)
    xgmax = last(degree_range(XiaoGimbutas(), Simplex{2}()))
    @test occursin("XiaoGimbutas", msg) && occursin("0–$xgmax", msg)
    @test occursin("0–∞", msg)
    e2 = try
        rule(XiaoGimbutas(), Simplex{2}(); degree = 60)
    catch err
        err
    end
    @test e2 isa NoRuleError
    @test occursin("Nearest: XiaoGimbutas degree $xgmax ($(npoints(XiaoGimbutas(), Simplex{2}(), xgmax)) points)",
                   sprint(showerror, e2))
    e3 = try
        rule(Simplex{2}(); degree = 9, positive = true, T = Rational{BigInt})
    catch err
        err
    end
    m3 = sprint(showerror, e3)
    @test occursin("has negative weights", m3)
    @test occursin("GrundmannMöller degree 1", m3)
    @test_throws NoRuleError rule(Simplex{2}(); degree = 5, family = XiaoGimbutas(), T = Rational{BigInt})
    @test_throws CR.NotYetImplemented Sphere{3}()
    @test_throws ArgumentError rule(Simplex{2}(); degree = -1)
    @test_throws ArgumentError rule(Simplex{2}(); npoints = 5)
end
