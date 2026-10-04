using CubatureRules, Test, LinearAlgebra
const CR = CubatureRules

# Rules for ∫_T f(y)/|y − x₀| dy on a triangle (notes/singular-bem.md). `check` verifies them
# against Dubiner polynomials integrated in polar coordinates about x₀; these tests add the
# closed forms, the cases that stress the construction, and the interface.

const TRI = Simplex((0, 0), (1, 0), (3//10, 8//10))

@testset "exact and sharp, with x₀ at a vertex, on an edge and inside" begin
    for (x0, npatch) in (((0, 0), 1), ((3//10, 8//10), 1), ((1//2, 0), 2), ((13//20, 4//10), 2),
                         ((1//3, 1//4), 3), ((0.41, 0.2), 3))
        dom = WeightedDomain(TRI, InverseDistance(x0))
        @test CR.npatches(dom) == npatch
        for d in (0, 1, 4, 9, 16)
            r = rule(dom; degree = d)
            m = max(1, cld(d + 1, 2))
            @test family(r) == "DuffyGauss"
            @test npoints(r) == npatch * m^2 == npoints(DuffyGauss(), dom, d)
            @test degree(r) == 2m - 1 >= d
            v = check(r)
            @test v.exact && v.sharp === true && v.positive && v.interior && v.weights_sum_ok
        end
    end
end

@testset "closed forms" begin
    # ∫ 1/|y| over the unit right triangle from its right-angle vertex: √2 log(1 + √2)
    r = rule(WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistance((0, 0))); degree = 9, digits = 50)
    @test abs(sum(weights(r)) - sqrt(big(2)) * log(1 + sqrt(big(2)))) < big(10.0)^-49
    # ∫ y₁/|y| over it: in polar coordinates ∫₀^{π/2} cos θ R(θ)²/2 dθ with R = 1/(cos θ + sin θ),
    # = (1/2) ∫₀^{π/2} cos θ / (1 + sin 2θ) dθ = 1/2 · (√2/2) log(1 + √2)
    @test abs(integrate(y -> y[1], r) - sqrt(big(2)) / 4 * log(1 + sqrt(big(2)))) < big(10.0)^-49
    # the measure of the weighted domain, as `check` uses it, against the rule
    dom = WeightedDomain(TRI, InverseDistance((1//3, 1//4)))
    @test abs(sum(weights(rule(dom; degree = 4, digits = 40))) - measure(dom)) < big(10.0)^-39
end

@testset "precision" begin
    dom = WeightedDomain(TRI, InverseDistance((1//2, 0)))
    for (kw, T) in (((T = Float32,), Float32), ((T = Float64,), Float64), ((digits = 50,), BigFloat))
        r = rule(dom; degree = 11, kw...)
        @test eltype(r) == T && check(r).exact
    end
    # the kernel's x₀ is read exactly: 0.3 is the binary number nearest 0.3
    @test CR.singular_patches(WeightedDomain(Simplex((0.0, 0.0), (0.6, 0.0), (0.0, 1.0)), InverseDistance((0.3, 0.0))))[2] === :edge
end

@testset "a thin sub-triangle: x₀ within 1e-6 of an edge" begin
    # the weight 1/√q(t) then has complex roots 1e-6 from [0, 1]; the sinh discretisation keeps
    # the construction cheap and the rule exact
    dom = WeightedDomain(TRI, InverseDistance((1//2, 1//1_000_000)))
    r = rule(dom; degree = 15, digits = 30)
    @test npoints(r) == 3 * 64
    v = check(r)
    @test v.exact && v.sharp === true && v.positive && v.interior
    # and a smooth integrand: the rule is not exact on e^{y₁}, but converges geometrically
    vals = [integrate(y -> exp(y[1]), rule(dom; degree = d, digits = 30)) for d in (11, 21, 31)]
    @test abs(vals[3] - vals[2]) < 1e-12 * abs(vals[3]) && abs(vals[3] - vals[2]) < abs(vals[2] - vals[1])
end

@testset "interface" begin
    dom = WeightedDomain(TRI, InverseDistance((1//3, 1//4)))
    @test CR.isreference(dom) && CR.reference(dom) === dom
    @test sprint(show, InverseDistance((1//3, 1//4))) == "InverseDistance((1//3, 1//4))"
    @test any(row -> row.family == "DuffyGauss", available(dom; degree = 6))
    r = rule(dom; degree = 6)
    @test occursin("inside", provenance(r).path[1]) && !isempty(provenance(r).citations)
    @test passed(verify(r))
    # refused: outside the triangle, a degenerate triangle, exact arithmetic
    @test_throws ArgumentError rule(WeightedDomain(TRI, InverseDistance((2, 0))); degree = 4)
    @test_throws ArgumentError rule(WeightedDomain(Simplex((0, 0), (1, 1), (2, 2)), InverseDistance((0, 0))); degree = 4)
    @test_throws NoRuleError rule(dom; degree = 4, T = Rational{BigInt})        # names DuffyGauss and why
    @test_throws ArgumentError InverseDistance((Inf, 0.0))
end
