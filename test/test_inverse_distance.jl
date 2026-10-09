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
    @test CR.kernel_geometry(WeightedDomain(Simplex((0.0, 0.0), (0.6, 0.0), (0.0, 1.0)), InverseDistance((0.3, 0.0)))).place === :edge
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

@testset "rounding is absorbed: a point meant to be on the triangle" begin
    # a third of the way along an edge, computed in Float64, is 7e-18 off it: moved onto it,
    # and the rule is built for that point
    T = Simplex((0.1, 0.0), (0.0, 0.3), (0.7, 0.9))
    dom = WeightedDomain(T, InverseDistance((0.1 + (0.0 - 0.1) / 3, 0.3 / 3)))
    (; place, moved) = CR.kernel_geometry(dom)
    @test place === :edge && 0 < moved < 1e-16
    r = rule(dom; degree = 7)
    @test npoints(r) == 2 * 16 && occursin("moved by", provenance(r).path[1])
    @test check(r).exact
    # exact input is never moved: a rational point 1e-30 from an edge is inside it
    (; place, moved) = CR.kernel_geometry(WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistance((1//2, 1//10^30))))
    @test place === :interior && moved == 0
end

@testset "triangles in space" begin
    # ordinary rules: the reference triangle's, mapped, with the area element
    T3 = SurfaceTriangle((0, 0, 0), (1, 0, 1), (0, 1, 1))
    @test measure(T3) ≈ sqrt(3) / 2
    for (kw, tol) in (((degree = 6,), 1e-14), ((degree = 6, digits = 40), 1e-39))
        r = rule(T3; kw...)
        v = check(r)
        @test v.exact && v.sharp === true && v.positive && v.interior
        @test abs(integrate(y -> y[3]^2, r) - sqrt(big(3)) / 4) < tol       # √3 ∫ (ξ₁ + ξ₂)² over the reference
    end
    # one reference rule for many triangles, mapped on the fly
    rref = rule(Simplex{2}(); degree = 6)
    f3(y) = y[3]^2 + y[1] * y[2]
    @test integrate(f3, rref, T3) ≈ integrate(f3, rule(T3; degree = 6)) rtol = 1e-14
    @test integrate(f3, rref, [T3, T3]) ≈ 2 * integrate(f3, rref, T3) rtol = 1e-14

    # the kernel: a planar triangle turned into space by a rational rotation has the planar
    # integrals, since the kernel depends only on distance
    Q = CR.SMatrix{3,3}(3//5, 0, -4//5, 0, 1, 0, 4//5, 0, 3//5) * CR.SMatrix{3,3}(1, 0, 0, 0, 3//5, 4//5, 0, -4//5, 3//5)
    shift = CR.SVector(1//2, -1//3, 2)
    lift3(p) = Q * CR.SVector(p[1], p[2], 0) + shift
    tri2 = [(0, 0), (1, 0), (3//10, 8//10)]
    f2(p) = 1 + p[1]^3 * p[2] - 2p[2]^4
    for x0 in ((1//3, 1//4), (1//2, 0), (0, 0))
        r2 = rule(WeightedDomain(Simplex(tri2...), InverseDistance(x0)); degree = 9, digits = 30)
        r3 = rule(WeightedDomain(SurfaceTriangle(lift3.(tri2)...), InverseDistance(lift3(x0))); degree = 9, digits = 30)
        @test npoints(r3) == npoints(r2)
        @test abs(integrate(y -> f2(transpose(Q) * (y - shift)), r3) - integrate(f2, r2)) < big(10.0)^-28
        v = check(r3)
        @test v.exact && v.sharp === true && v.positive && v.interior
    end

    # a centroid computed in Float64 is off the plane by rounding: projected onto it
    v = ((0.1, 0.2, 0.3), (1.3, -0.4, 0.7), (0.2, 1.1, -0.5))
    c = Tuple((CR.SVector(v[1]) + CR.SVector(v[2]) + CR.SVector(v[3])) / 3)
    r = rule(WeightedDomain(SurfaceTriangle(v...), InverseDistance(c)); degree = 9)
    @test npoints(r) == 3 * 25 && check(r).exact
    @test_throws ArgumentError SurfaceTriangle((0, 0, 0), (1, 1, 1), (2, 2, 2))
    @test_throws ArgumentError SurfaceTriangle((0, 0), (1, 0), (0, 1))
end

@testset "near-singular: x₀ off the triangle" begin
    # above a triangle in space, from far to very close: exact, and tending to the singular
    # integral at the foot as x₀ comes down, the gap being about 2πH f(foot)
    T = SurfaceTriangle((0, 0, 0), (1, 0, 0), (3//10, 8//10, 0))
    f(y) = 1 + y[1]^2 * y[2] - y[2]^3
    Is = integrate(f, rule(WeightedDomain(T, InverseDistance((1//3, 1//4, 0))); degree = 9, digits = 30))
    for H in (1//10^4, 1//10^8)
        r = rule(WeightedDomain(T, InverseDistance((1//3, 1//4, H))); degree = 9)
        @test family(r) == "DuffySinh"
        @test abs(integrate(f, r) - Is + 2π * H * f((1//3, 1//4))) < 50H^2 + 1e-14
        # `check` costs seconds here (its reference integrals are refined to 128 bits): once
        if H == 1//10^8
            v = check(r)
            @test v.exact && v.sharp === true && v.positive && v.interior
        end
    end
    # the point count does not grow as x₀ comes closer: each radial line has the Gauss rule of
    # its own weight (and further away it falls, the angular substitution scaled to the height)
    n12, n6, n1 = (npoints(rule(WeightedDomain(T, InverseDistance((1//3, 1//4, H))); degree = 9))
                   for H in (1//10^12, 1//10^6, 1//10))
    @test n12 <= n6 && n1 <= n6
    # in the plane, just outside an edge
    r = rule(WeightedDomain(TRI, InverseDistance((1//2, -1//1000))); degree = 7)
    v = check(r)
    @test v.exact && v.sharp === true && v.positive && v.interior
    # at 30 digits, against the 2πH asymptotics
    r = rule(WeightedDomain(T, InverseDistance((1//3, 1//4, 1//10^12))); degree = 9, digits = 30)
    @test abs(integrate(f, r) - Is + 2π * big(1//10^12) * f((1//3, 1//4))) < big(10.0)^-21
    # far away, an ordinary rule of high degree agrees
    far = WeightedDomain(TRI, InverseDistance((2, 2)))
    @test integrate(f, rule(far; degree = 7)) ≈
          integrate(y -> f(y) / sqrt((y[1] - 2)^2 + (y[2] - 2)^2), rule(TRI; degree = 40)) rtol = 1e-14
end

@testset "interface" begin
    dom = WeightedDomain(TRI, InverseDistance((1//3, 1//4)))
    @test CR.isreference(dom) && CR.reference(dom) === dom
    @test sprint(show, InverseDistance((1//3, 1//4))) == "InverseDistance((1//3, 1//4))"
    @test any(row -> row.family == "DuffyGauss", available(dom; degree = 6))
    r = rule(dom; degree = 6)
    @test occursin("inside", provenance(r).path[1]) && !isempty(provenance(r).citations)
    @test passed(verify(r))
    # a point off the triangle is a near-singular one, for DuffySinh
    @test family(rule(WeightedDomain(TRI, InverseDistance((2, 0))); degree = 4)) == "DuffySinh"
    # refused: a degenerate triangle, exact arithmetic
    @test_throws ArgumentError rule(WeightedDomain(Simplex((0, 0), (1, 1), (2, 2)), InverseDistance((0, 0))); degree = 4)
    @test_throws NoRuleError rule(dom; degree = 4, T = Rational{BigInt})        # names DuffyGauss and why
    @test_throws ArgumentError InverseDistance((Inf, 0.0))
end
