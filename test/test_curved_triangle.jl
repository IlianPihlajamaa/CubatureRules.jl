using CubatureRules, LinearAlgebra, StaticArrays, Test
const CR = CubatureRules

# The weakly singular kernel on a curved (quadratic) triangle, as the weight J(ξ)/|χ(ξ) − χ(ξ₀)|
# on the reference triangle (DuffyCurved). Besides `check`, which integrates the Dubiner
# polynomials in polar coordinates about ξ₀, a straight-sided element against DuffyGauss on the
# same flat triangle, and the geometry against finite differences.

s = 1 / sqrt(2)
Γ = QuadraticTriangle((1, 0, 0), (0, 1, 0), (0, 0, 1), (s, s, 0), (0, s, s), (s, 0, s))   # an octant of the sphere

@testset "geometry" begin
    @test Γ((0, 0)) ≈ [1, 0, 0] && Γ((1 // 2, 0)) ≈ [s, s, 0] && Γ((0, 1 // 2)) ≈ [s, 0, s]
    ξ, h = SVector(0.2, 0.3), 1e-6
    A = CR.jacobian_matrix(Γ, ξ)
    @test A[:, 1] ≈ (Γ(ξ + SVector(h, 0)) - Γ(ξ - SVector(h, 0))) / 2h rtol = 1e-8
    @test A[:, 2] ≈ (Γ(ξ + SVector(0, h)) - Γ(ξ - SVector(0, h))) / 2h rtol = 1e-8
    # χ(ξ₀ + s v) − χ(ξ₀) = s Dχ(ξ₀) v + s² Q(v), exactly for a quadratic map
    v, σ = SVector(0.3, -0.1), 0.7
    @test Γ(ξ + σ * v) - Γ(ξ) ≈ σ * A * v + σ^2 * CR.second_term(Γ, v) rtol = 1e-13
    # a straight-sided element has the area element of its flat triangle
    P1, P2, P3 = (0.0, 0.0, 0.0), (1.0, 0.2, 0.1), (0.3, 0.9, -0.2)
    mid(a, b) = (a .+ b) ./ 2
    Γf = QuadraticTriangle(P1, P2, P3, mid(P1, P2), mid(P2, P3), mid(P3, P1))
    @test CR.area_element(Γf, ξ) ≈ norm(cross(collect(P2 .- P1), collect(P3 .- P1)))
end

@testset "rules: exact, sharp, positive, interior" begin
    for (ξ0, d) in (((0, 0), 5), ((1 // 2, 0), 3), ((1 // 5, 3 // 10), 5))
        dom = WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ, ξ0))
        r = rule(dom; degree = d)
        v = check(r)
        @test family(r) == "DuffyCurved" && degree(r) >= d
        @test passed(v) && v.exact && v.sharp === true && v.positive && v.interior
    end
    dom = WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ, (1 // 5, 3 // 10)))
    r = rule(dom; degree = 4, digits = 30)
    @test abs(sum(weights(r)) - setprecision(() -> CR.measure(dom), BigFloat, 256)) < 1e-28
end

@testset "a straight-sided element is a flat triangle" begin
    P1, P2, P3 = (0.0, 0.0, 0.0), (1.0, 0.2, 0.1), (0.3, 0.9, -0.2)
    mid(a, b) = (a .+ b) ./ 2
    Γf = QuadraticTriangle(P1, P2, P3, mid(P1, P2), mid(P2, P3), mid(P3, P1))
    f(y) = 1 + y[1] * y[2] - y[3]^2
    for ξ0 in ((0, 0), (1 // 2, 1 // 2), (1 // 4, 1 // 3))
        x0 = Γf(ξ0)
        rc = rule(WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γf, ξ0)); degree = 8)
        rf = rule(WeightedDomain(SurfaceTriangle(P1, P2, P3), InverseDistance(Tuple(x0))); degree = 8)
        @test integrate(ξ -> f(Γf(ξ)), rc) ≈ integrate(f, rf) rtol = 1e-13
    end
end
