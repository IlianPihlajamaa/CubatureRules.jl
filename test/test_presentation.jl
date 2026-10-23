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
    # a rule of the package's own table: the package first, then the class of rules it follows
    @test startswith(b, "@misc{Pihlajamaa_CubatureRules,")
    @test occursin("@article{XiaoGimbutas2010,", b) && occursin("doi = {10.1016/j.camwa.2009.10.027}", b)
    @test !occursin("Compared with", b)
    @test occursin("Xiao, H., & Gimbutas, Z. (2010).", cite(r; style = :apa))
    @test occursin("doi:10.1016", cite(r; style = :plain))
    @test_throws ArgumentError cite(r; style = :mla)
    @test occursin("@book{Stroud1971", cite(rule(ConicalProduct(), Simplex{2}(); degree = 3)))
    # published rules a table was only compared with are no sources: listed apart, after them
    c = rule(FullySymmetric(), Orthotope{3}(); degree = 9)
    ckeys(cs) = [x.key for x in cs]
    @test ckeys(provenance(c).citations) == ["Pihlajamaa_CubatureRules"]
    @test ckeys(provenance(c).comparisons) == ["WitherdenVincent2015"]
    bc = cite(c)
    head, wv = findfirst("% Compared with, not sources of this rule:", bc), findfirst("WitherdenVincent2015", bc)
    @test head !== nothing && wv !== nothing && first(wv) > last(head)
    @test !occursin("WitherdenVincent2015", cite(c; comparisons = false))
    @test occursin("Compared with, not sources", cite(c; style = :plain))
    shown = sprint(show, MIME"text/plain"(), c)
    @test occursin("reference : CubatureRules.jl (2026)", shown) && occursin("compared  : Witherden & Vincent (2015)", shown)
    @test occursin("compared  : WitherdenVincent2015", sprint(show, MIME"text/plain"(), provenance(c)))
    # Witherden & Vincent have no 4-cube rules; the sphere's degree 19-131 grew from the
    # package's own degree-17 rule, and its 133-201 from Lebedev & Laikov's degree-131 one
    @test isempty(provenance(rule(FullySymmetric(), Orthotope{4}(); degree = 5)).comparisons)
    s41 = provenance(rule(LebedevRule(), Sphere{3}(); degree = 41))
    @test ckeys(s41.citations) == ["Pihlajamaa_CubatureRules", "Lebedev1976"] && ckeys(s41.comparisons) == ["LebedevLaikov1999"]
    s133 = provenance(rule(LebedevRule(), Sphere{3}(); degree = 133))
    @test ckeys(s133.citations) == ["Pihlajamaa_CubatureRules", "Lebedev1976", "LebedevLaikov1999"] && isempty(s133.comparisons)
    # a product rule carries its parts' comparisons along with their citations
    b41 = provenance(rule(BallProduct(LebedevRule()), Ball{3}(); degree = 41))
    @test "Lebedev1976" in ckeys(b41.citations) && ckeys(b41.comparisons) == ["LebedevLaikov1999"]
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
