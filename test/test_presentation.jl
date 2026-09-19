using CubatureRules, Test
const CR = CubatureRules

@testset "show" begin
    r = rule(Simplex{2}(); degree = 17, digits = 40)
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("QuadratureRule{2,BigFloat} on Simplex{2}()", s)
    @test occursin("family    : XiaoGimbutas (seeded, Newton-refined)", s)
    @test occursin("polynomial degree 17", s)
    @test occursin("60, all interior", s)
    @test occursin("all positive", s)
    @test occursin("est. cond(J)", s)
    @test occursin("residual  :", s) && occursin("defining equations", s)   # names which residual
    @test occursin("Xiao & Gimbutas (2010)", s)
    @test occursin("XiaoGimbutas, 60 points", sprint(show, r))
    g = rule(GrundmannMöller(), Simplex{2}(); degree = 5, T = Rational{BigInt})
    sg = sprint(show, MIME"text/plain"(), g)
    @test occursin("exact rational", sg) && occursin("some negative", sg)
    @test occursin("PASSED", sprint(show, MIME"text/plain"(), check(g)))
    @test occursin("Certificate", sprint(show, MIME"text/plain"(), certificate(r)))
    @test occursin("seed", sprint(show, MIME"text/plain"(), provenance(r)))
    @test occursin("Candidates", sprint(show, MIME"text/plain"(), available(Simplex{2}(); degree = 5)))
end

@testset "cite" begin
    r = rule(Simplex{2}(); degree = 10)
    b = cite(r)
    @test startswith(b, "@article{XiaoGimbutas2010,")
    @test occursin("doi = {10.1016/j.camwa.2009.10.027}", b)
    @test occursin("Xiao, H., & Gimbutas, Z. (2010).", cite(r; style = :apa))
    @test occursin("doi:10.1016", cite(r; style = :plain))
    @test_throws ArgumentError cite(r; style = :mla)
    @test occursin("@book{Stroud1971", cite(rule(ConicalProduct(), Simplex{2}(); degree = 3)))
end

@testset "construction benchmark" begin
    rows = benchmark_construction(Simplex{2}(); degrees = [3, 60])
    @test any(r -> r.family == "XiaoGimbutas" && r.degree == 3 && r.status === :completed, rows)
    @test any(r -> r.family == "XiaoGimbutas" && r.degree == 60 && r.status === :unsupported, rows)
    @test all(r -> r.status !== :completed || r.seconds >= 0, rows)
    rq = benchmark_construction(Simplex{2}(); degrees = [3], T = Rational{BigInt})
    @test any(r -> r.family == "ConicalProduct(GaussJacobi)" && r.status === :unsupported, rq)
end
