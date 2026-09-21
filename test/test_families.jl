using CubatureRules, Test, StaticArrays
const CR = CubatureRules

@testset "Gauss–Jacobi" begin
    for n in (1, 2, 5, 20, 64)
        r = rule(GaussLegendre(), Interval(); npoints = n)
        @test npoints(r) == n
        @test degree(r) == 2n - 1
        @test issorted(nodes(r))
        @test sum(weights(r)) ≈ 2
        @test passed(check(r))
    end
    # known values
    r3 = rule(GaussLegendre(), Interval(); npoints = 3, digits = 50)
    @test abs(nodes(r3)[3] - sqrt(big(3) / 5)) < big(10.0)^-49
    @test abs(weights(r3)[2] - big(8) / 9) < big(10.0)^-49
    # Jacobi weights live on a weighted domain
    wd = WeightedDomain(Interval(), JacobiWeight(2, 1))
    rj = rule(wd; degree = 9)
    @test domain(rj) == wd
    @test sum(weights(rj)) ≈ Float64(measure(wd))
    @test passed(check(rj))
    rj2 = rule(WeightedDomain(Interval(), JacobiWeight(-0.5, -0.5)); degree = 11, digits = 40)
    @test abs(sum(weights(rj2)) - big(pi)) < big(10.0)^-38          # Chebyshev: ∫ = π
    @test passed(check(rj2))
    @test_throws ArgumentError GaussJacobi(-1, 0)
end

@testset "Grundmann–Möller, exact rationals, any dimension" begin
    for D in 1:4, d in (1, 3, 5, 7)
        D == 1 && continue
        r = rule(GrundmannMöller(), Simplex{D}(); degree = d, T = Rational{BigInt})
        @test eltype(r) == Rational{BigInt}
        @test sum(weights(r)) == 1 // factorial(D)
        @test npoints(r) == npoints(GrundmannMöller(), Simplex{D}(), d)
        v = check(r)
        @test v.exact
        @test v.sharp === true
        @test v.symmetric === true
        @test v.positive == (d == 1)
    end
    # even degree requests round up to the next odd degree
    @test degree(rule(GrundmannMöller(), Simplex{2}(); degree = 4)) == 5
    @test GrundmannMoeller === GrundmannMöller
    # floating output is the exact rule rounded once
    rf = rule(GrundmannMöller(), Simplex{3}(); degree = 5, digits = 40)
    re = rule(GrundmannMöller(), Simplex{3}(); degree = 5, T = Rational{BigInt})
    @test all(abs.(weights(rf) .- weights(re)) .< big(10.0)^-40)
end

@testset "conical product, any dimension" begin
    for D in (2, 3, 4), d in (1, 4, 9)
        r = rule(ConicalProduct(), Simplex{D}(); degree = d)
        @test npoints(r) == cld(d + 1, 2)^D
        @test sum(weights(r)) ≈ 1 / factorial(D)
        v = check(r)
        @test v.exact && v.sharp === true && v.positive && v.interior
    end
    @test_throws NoRuleError rule(ConicalProduct(), Simplex{2}(); degree = 3, T = Rational{BigInt})
end

@testset "Xiao–Gimbutas: every shipped degree verifies (Float64)" begin
    # counts published by Xiao & Gimbutas (2010), Table 1, column n6. A shipped rule may
    # differ either way — node elimination has found smaller ones — but never by much, and
    # any difference must be recorded in the entry's note.
    counts = [1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79,
              87, 96, 103, 112, 120, 130, 141, 150, 159, 171, 181, 193, 204, 214, 228,
              243, 252, 267, 282, 295, 309, 324, 339, 354, 370, 385, 399, 423, 435, 453]
    dmax = last(degree_range(XiaoGimbutas(), Simplex{2}()))
    @test dmax >= 20
    for d in 1:dmax
        r = rule(XiaoGimbutas(), Simplex{2}(); degree = d)
        e = CR.xg_entry_for(d)
        @test npoints(r) <= counts[e.degree] + 3
        @test npoints(r) == counts[e.degree] || !isempty(e.note)
        v = check(r)
        @test v.exact
        @test v.positive && v.interior && v.symmetric === true
        @test v.sharp !== false
        @test certificate(r).residual < 1e-14
    end
end

@testset "Xiao–Gimbutas at arbitrary precision (the cornerstone)" begin
    r = rule(Simplex{2}(); degree = 20, digits = 200)
    @test family(r) == "XiaoGimbutas"
    @test npoints(r) == 79
    @test precision(first(weights(r))) == CR.digits_to_bits(200)
    v = check(r)
    @test v.exact                       # exact at degree 20 …
    @test v.sharp === true              # … and not at 21
    @test v.max_residual < big(10.0)^-195
    @test v.symmetric === true && v.positive && v.interior
    c = certificate(r)
    @test c.residual < big(10.0)^-199
    @test c.iterations >= 1
    # refinement agrees with the Float64 rule to Float64 precision
    r64 = rule(Simplex{2}(); degree = 20)
    @test maximum(abs.(Float64.(weights(r)) .- weights(r64))) < 1e-16
    # higher precision is consistent with lower precision
    r100 = rule(Simplex{2}(); degree = 20, digits = 100)
    @test maximum(abs.(weights(r) .- weights(r100))) < big(10.0)^-99
end

@testset "FullySymmetric tetrahedra: every shipped degree verifies" begin
    rg = degree_range(FullySymmetric(), Simplex{3}())
    @test last(rg) >= 8
    for d in 1:last(rg)
        r = rule(FullySymmetric(), Simplex{3}(); degree = d)
        @test degree(r) >= d
        @test sum(weights(r)) ≈ 1 / 6
        v = check(r)
        @test v.exact
        @test v.positive && v.interior && v.symmetric === true
        @test v.sharp !== false
        @test certificate(r).residual < 1e-14
    end
    # arbitrary precision, and the selector prefers it on the tetrahedron
    r = rule(Simplex{3}(); degree = 8, digits = 60)
    @test family(r) == "FullySymmetric"
    v = check(r)
    @test v.exact && v.sharp === true && v.max_residual < big(10.0)^-55
    @test isempty(candidates(FullySymmetric, Simplex{2}(), PolynomialDegree(3)))
end

@testset "shipped symmetric rules have distinct nodes" begin
    # an orbit collapsing onto another would make the point count meaningless
    for (fam, dom) in ((XiaoGimbutas(), Simplex{2}()), (FullySymmetric(), Simplex{3}()))
        for d in 1:last(degree_range(fam, dom))
            x = nodes(rule(fam, dom; degree = d))
            sep = minimum(maximum(abs, x[i] - x[j]) for i in eachindex(x) for j in (i + 1):length(x); init = Inf)
            @test sep > 1e-6
        end
    end
end

@testset "Newton–Cotes, exact rationals" begin
    # the classical small rules, on [-1, 1]
    @test CR.newton_cotes_exact(2, false)[2] == [1, 1]                       # trapezoid
    @test CR.newton_cotes_exact(3, false)[2] == [1 // 3, 4 // 3, 1 // 3]     # Simpson
    @test CR.newton_cotes_exact(5, false)[2] == [7 // 45, 32 // 45, 4 // 15, 32 // 45, 7 // 45]  # Boole
    @test CR.newton_cotes_exact(1, true) == ([0 // 1], [2 // 1])             # midpoint
    for open in (false, true), d in (0, 1, 2, 3, 5, 7, 9)
        f = NewtonCotes(open)
        r = rule(f, Interval(); degree = d, T = Rational{BigInt})
        v = check(r)
        @test degree(r) >= d
        @test sum(weights(r)) == 2
        @test v.exact && v.sharp === true && v.symmetric === true
        @test v.interior == open
        @test v.positive == properties(f, Interval(), d).positive
        @test eltype(r) == Rational{BigInt}
    end
    # floating output is the exact rule rounded once
    rf = rule(NewtonCotes(), Interval(); degree = 7, digits = 50)
    re = rule(NewtonCotes(), Interval(); degree = 7, T = Rational{BigInt})
    @test all(abs.(weights(rf) .- weights(re)) .< big(10.0)^-49)
    @test passed(check(rf))
    # Gauss needs fewer points, so the selector prefers it
    @test family(rule(Interval(); degree = 7)) == "GaussJacobi"
    @test npoints(NewtonCotes(), Interval(), 7) == 7 > npoints(GaussLegendre(), Interval(), 7)
end

@testset "Fejér" begin
    for kind in (1, 2), d in (0, 1, 2, 5, 8, 15)
        f = Fejer(kind)
        r = rule(f, Interval(); degree = d)
        v = check(r)
        @test degree(r) >= d
        @test npoints(r) == (max(1, isodd(d) ? d : d + 1))
        @test sum(weights(r)) ≈ 2
        @test v.exact && v.sharp === true
        @test v.positive && v.interior && v.symmetric === true     # positive at every order
    end
    # arbitrary precision from the closed form
    r = rule(Fejer(2), Interval(); degree = 21, digits = 60)
    @test passed(check(r))
    @test abs(sum(weights(r)) - 2) < big(10.0)^-58
    @test abs(integrate(exp, r) - (exp(big(1)) - exp(big(-1)))) < big(10.0)^-25   # not exact: exp is not a polynomial
    @test_throws ArgumentError Fejer(3)
    @test_throws NoRuleError rule(Fejer(1), Interval(); degree = 5, T = Rational{BigInt})
    # Gauss still wins on node count, so the selector prefers it
    @test family(rule(Interval(); degree = 9)) == "GaussJacobi"
end
