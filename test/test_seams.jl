using CubatureRules, LinearAlgebra, Random, StaticArrays, Test
import CubatureRules: ⊗   # public, not exported
const CR = CubatureRules

# Where components meet: what a request on one domain gets from a family made for another,
# what a domain's number type does to a rule, what a rule passed between two operations keeps.

@testset "an explicit family answers for its own domain only" begin
    # GaussJacobi(1, 1) integrates (1 − x)(1 + x) f; asked for the plain interval, it must say so
    @test_throws NoRuleError rule(Interval(); degree = 5, family = GaussJacobi(1, 1))
    @test_throws NoRuleError rule(GaussLaguerre(1.0), LaguerreRay(0); degree = 5)
    @test sum(weights(rule(GaussLaguerre(1.0), LaguerreRay(1.0); degree = 5))) ≈ 1
    @test sum(weights(rule(GaussJacobi(1, 1), WeightedDomain(Interval(), JacobiWeight(1, 1)); degree = 5))) ≈ 4 / 3
end

@testset "the domain's coordinates do not set the rule's number type" begin
    r = rule(Interval(big(0), big(1)); degree = 5)
    @test numbertype(r) === Float64 && passed(check(r))
    t = Simplex((big(0), big(0)), (big(1), big(0)), (big(0), big(1)) ./ 3)
    rt = rule(t; degree = 4)
    @test numbertype(rt) === Float64 && passed(check(rt))
    rb = rule(t; degree = 4, digits = 30)
    @test numbertype(rb) === BigFloat && precision(first(weights(rb))) == CR.digits_to_bits(30) && passed(check(rb))
    # a Float32 request on a Float64 mesh is a Float32 rule
    @test numbertype(rule(Interval(0.0, 2.0); degree = 5, T = Float32)) === Float32
end

@testset "domains of mixed coordinates, and degenerate ones" begin
    @test Simplex((0, 0), (1, 0), (0.5, 1)) isa Simplex{2,Float64}
    @test Orthotope((0, 0.5), (1, 2)) isa Orthotope{2,Float64}
    @test Ball((0, 0.5, 1), 2) isa Ball{3,Float64}
    @test_throws ArgumentError Simplex((0, 0), (1, 1), (2, 2))
    @test_throws ArgumentError Simplex((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0))
    # exact vertices, exact measure
    m = measure(Simplex((0, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1)))
    @test m isa Rational && m == 1 // 3
    @test measure(Simplex((0, 0), (1 // 2, 0), (0, 1))) == 1 // 4
end

@testset "weights that live on one kind of domain" begin
    @test_throws ArgumentError WeightedDomain(Simplex{2}(), JacobiWeight(1, 0))
    @test_throws ArgumentError WeightedDomain(Interval(), InverseDistance((0.0, 0.0)))
    @test_throws ArgumentError WeightedDomain(Interval(), CR.GaussianWeight())
    # a Jacobi weight on [a, b] is (b − x)^α (x − a)^β
    d = WeightedDomain(Interval(0, 2), JacobiWeight(1, 0))
    @test measure(d) == 2
    r = rule(d; degree = 5)
    @test passed(check(r))
    @test integrate(x -> x^2, r) ≈ 4 / 3                    # ∫₀² x² (2 − x) dx
    @test integrate(x -> x^2, rule(WeightedDomain(Interval(0, 1), JacobiWeight(0, 2)); degree = 7)) ≈ 1 / 5
end

@testset "embedded Gauss–Kronrod on any interval" begin
    e = embedded(GaussKronrod(), Interval(0, 1); degree = 11)
    res = integrate(exp, e; error = true)
    @test res.value ≈ exp(1) - 1 && 0 < res.error_estimate < 1e-5
    @test sum(e.coarse_weights) ≈ 1
end

@testset "products and transforms keep the rules' precision" begin
    a = rule(Interval(); degree = 9, digits = 100)
    p = setprecision(() -> a ⊗ a, BigFloat, 64)
    @test precision(first(weights(p))) == precision(first(weights(a)))
    @test passed(check(p))
    t = setprecision(() -> CR.transform(a, x -> x / 2, x -> 1 // 2; domain = Interval(-1 // 2, 1 // 2)), BigFloat, 64)
    @test precision(first(weights(t))) == precision(first(weights(a)))
end

@testset "what show and check say about the weights" begin
    r = rule(Simplex((0, 0), (1 // 1000, 0), (0, 1 // 1000)); degree = 4, digits = 40)
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("Σw = 5e-07", s)
    # a weight that underflows to zero is not negative
    l = rule(LaguerreRay(0); degree = 399)
    @test any(iszero, weights(l))
    @test check(l).positive
    @test occursin("non-negative, 1 zero after rounding to Float64", sprint(show, MIME"text/plain"(), l))
end

@testset "rule from several tasks at once" begin
    CR.selection_warnings!(true)
    rs = fetch.([Threads.@spawn rule(Simplex{2}(); degree = d) for d in 1:12])
    @test all(passed ∘ check, rs)
    @test all(i -> npoints(rs[i]) == npoints(rule(Simplex{2}(); degree = i)), 1:12)
end

@testset "a rule is a collection of (node, weight) pairs" begin
    r = rule(Interval(); degree = 5)
    @test length(r) == npoints(r) == length(collect(r))
    @test sum(w * x^2 for (x, w) in r) ≈ 2 / 3
    @test eltype(r) == Tuple{Float64,Float64} && numbertype(r) === Float64
end

@testset "requests checked where they are made" begin
    @test numbertype(rule(Interval(); degree = 5, digits = 30.0)) === BigFloat
    @test_throws ArgumentError rule(Interval(); degree = 5, digits = 30.5)
    @test occursin("Rational{BigInt}", sprint(showerror, try rule(Interval(); degree = 5, T = Int) catch e e end))
    # npoints is the fewest a rule can have when its construction decides the size
    k = WeightedDomain(Interval(), FinitePart(9 // 10))
    @test !CR.exact_npoints(CR.SingularGauss(), k, 8)
    r = rule(CR.SingularGauss(), k; npoints = 5)
    @test npoints(r) == 5 || occursin("the construction needed $(npoints(r))", provenance(r).selection)
end

@testset "families and their claims" begin
    # Gauss–Jacobi with α = β is symmetric, and says so, and check confirms it
    r = rule(GaussJacobi(1 // 2, 1 // 2), WeightedDomain(Interval(), JacobiWeight(1 // 2, 1 // 2)); degree = 9, digits = 30)
    @test provenance(r).symmetry isa CR.Reflection && check(r).symmetric === true
    @test CR.properties(GaussJacobi(1, 0), Interval(), 5).symmetry isa CR.NoSymmetry
    @test_throws ArgumentError subdivide(rule(WeightedDomain(Interval(), JacobiWeight(1, 0)); degree = 3), 2)
end

@testset "citations" begin
    q = CR.cite(rule(GaussKronrod(), Interval(); degree = 11))
    @test occursin("@misc{Johnson_QuadGK", q) && occursin("url = {https://github.com/JuliaMath/QuadGK.jl}", q)
    @test !occursin("@book{Johnson_QuadGK", q)
    @test CR._bibescape("A & B_1") == "A \\& B\\_1"
    r = CR.transform(rule(Interval(); degree = 3), x -> x, x -> 1)
    t = CR.QuadratureRule(nodes(r), weights(r), domain(r), exactness(r),
                          CR.Provenance(provenance(r); citations = CR.Citation[]), nothing)
    @test occursin("no references", cite(t))
end

@testset "orthonormal verification bases" begin
    for (dom, d, name) in ((Disk(), 9, "Zernike"), (Ball{3}(), 9, "ball polynomials"),
                           (GaussianSpace(3), 7, "Hermite"), (Simplex{4}(), 5, "4-simplex"))
        r = rule(dom; degree = d, digits = 30)
        v = check(r)
        @test passed(v) && occursin(name, v.basis)
    end
    # a rule that is wrong in one weight fails
    r = rule(FullySymmetric(), Ball{3}(); degree = 9, digits = 30)
    w = copy(weights(r))
    w[3] *= 1 + big"1e-20"
    @test !passed(check(CR.QuadratureRule(nodes(r), w, domain(r), exactness(r), provenance(r), certificate(r))))
    # exact rules on the 4-simplex are checked exactly
    v = check(rule(GrundmannMöller(), Simplex{4}(); degree = 5, T = Rational{BigInt}))
    @test passed(v) && v.max_residual == 0
end

@testset "rule_hash does not read show" begin
    r = rule(Simplex((0, 0), (1, 0), (0, 1 // 2)); degree = 3)
    io = IOBuffer()
    CR._canonical(io, domain(r))
    s = String(take!(io))
    @test occursin("Simplex(", s) && occursin("0.5", s) && !occursin(sprint(show, domain(r)), s)
end

@testset "the condition number behind the guard bits" begin
    # no LAPACK, so the same on every platform; against the SVD to within a few percent
    rng = Random.Xoshiro(1)
    for (m, n, s) in ((3, 3, 0), (8, 5, 4), (20, 20, 10), (20, 8, 6))
        A = randn(rng, m, n) * Diagonal(10.0 .^ range(0, -s; length = n))
        σ = svdvals(A)
        @test CR.estimate_cond(A) ≈ σ[1] / σ[end] rtol = 0.05
        @test CR.estimate_cond(permutedims(A)) ≈ CR.estimate_cond(A) rtol = 0.05
    end
end
