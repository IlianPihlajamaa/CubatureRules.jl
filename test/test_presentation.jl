using CubatureRules, Test
import CubatureRules: Certificate, benchmark_construction   # public, not exported
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
    beyond = last(CubatureRules.degree_range(XiaoGimbutas(), Simplex{2}())) + 1   # past the seeded table
    rows = benchmark_construction(Simplex{2}(); degrees = [3, beyond])
    @test any(r -> r.family == "XiaoGimbutas" && r.degree == 3 && r.status === :completed, rows)
    @test any(r -> r.family == "XiaoGimbutas" && r.degree == beyond && r.status === :unsupported, rows)
    @test all(r -> r.status !== :completed || r.seconds >= 0, rows)
    rq = benchmark_construction(Simplex{2}(); degrees = [3], T = Rational{BigInt})
    @test any(r -> r.family == "ConicalProduct(GaussJacobi)" && r.status === :unsupported, rq)
end

@testset "error messages are wrapped, but not inside code" begin
    s = "a long explanation, " ^ 8 * "then `rule(ExpSinh(4), HalfLine())` and Interval(0, 1), and [0, 1]."
    m = CubatureRules.wrap_message(s; width = 40)
    lines = split(m, '\n')
    @test textwidth(lines[1]) <= 33 && all(l -> textwidth(l) <= 40, lines)
    @test any(l -> occursin("`rule(ExpSinh(4), HalfLine())`", l), lines)       # code is never split
    @test any(l -> occursin("Interval(0, 1),", l), lines) && any(l -> occursin("[0, 1].", l), lines)
    @test join(lines, ' ') == s                                                # only spaces became breaks
    # indented lines keep their indentation, and short lines are left alone
    @test CubatureRules.wrap_message("short\n  " * "item " ^ 20; width = 30) ==
          "short\n" * join(["  item item item item item", "    item item item item item",
                            "    item item item item item", "    item item item item item"], '\n')
    # what the package's errors print
    e = try rule(Simplex{2}(); degree = 9, positive = true, T = Rational{BigInt}) catch err err end
    shown = split(sprint(showerror, e), '\n')
    @test textwidth(shown[1]) <= 85 && all(l -> textwidth(l) <= 92, shown)
end
