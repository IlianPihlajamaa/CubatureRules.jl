using CubatureRules, LinearAlgebra, StaticArrays, Test
const CR = CubatureRules

# The weakly singular kernel on a curved (quadratic) triangle, as the weight J(ξ)/|χ(ξ) − χ(ξ₀)|
# on the reference triangle (DuffyCurved). Besides `check`, which integrates the Dubiner
# polynomials in polar coordinates about ξ₀, a straight-sided element against DuffyGauss on the
# same flat triangle, and the geometry against finite differences. With x₀ off the element, the
# nearest point against a grid, the limit onto the surface, and DuffySinh on a flat element.

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

@testset "near-singular: x₀ off the element" begin
    # x₀ given in space: cut at its nearest point χ(ξc), each ray's weight s J/|χ − x₀| from a
    # discretisation scaled to the distance
    ξ = (1 // 5, 3 // 10)
    y = Γ(ξ)
    A = CR.jacobian_matrix(Γ, ξ)
    n = normalize(cross(A[:, 1], A[:, 2]))
    grid_distance(x0) = minimum(norm(Γ((i / 300, j / 300)) - collect(x0)) for i in 0:300 for j in 0:(300 - i))
    # above and below (inside the sphere) the element, beyond an edge; inside the sphere near an
    # edge, where the nearest point is inside the element; beyond a vertex. `check` (polar
    # integrals about ξc, slow for a point this close) on two, the weight sum on one
    for (x0, ξc, verify) in ((y + 1e-2 * n, (0.2, 0.3), :check), (y - 1e-2 * n, (0.2, 0.3), :none),
                             ((0.75, 0.75, -0.05), (0.5, 0.0), :check), ((0.6, 0.6, -0.01), (0.48781, 0.02438), :none),
                             ((1.01, -0.01, -0.01), (0.0, 0.0), :measure))
        dom = WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ; x0))
        g = CR.curved_geometry(dom)
        @test g.place === :near && norm(Float64.(g.x0) - collect(ξc)) < 1e-5
        @test g.H <= grid_distance(x0) + 1e-12
        r = rule(dom; degree = 3)
        @test family(r) == "DuffyCurved" && degree(r) >= 3 && all(>(0), weights(r))
        if verify === :check
            v = check(r)
            @test passed(v) && v.exact && v.sharp === true && v.interior
        elseif verify === :measure
            @test sum(weights(r)) ≈ setprecision(() -> CR.measure(dom), BigFloat, 128) rtol = 1e-14
        end
    end
    # as x₀ comes down onto the surface the integral tends to the singular one, linearly
    f(ξ) = 1 + ξ[1] - ξ[2]^2
    I0 = integrate(f, rule(WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ, ξ)); degree = 3))
    dI = [abs(integrate(f, rule(WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ; x0 = y + H * n)); degree = 3)) - I0)
          for H in (1e-4, 1e-8)]
    @test 0.9e4 < dI[1] / dI[2] < 1.1e4
    # a point given in space on the surface, to rounding, is the singular case
    @test CR.curved_geometry(WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ; x0 = y))).place === :interior
    # a straight-sided element: against DuffySinh on the flat triangle
    P1, P2, P3 = (0.0, 0.0, 0.0), (1.0, 0.2, 0.1), (0.3, 0.9, -0.2)
    mid(a, b) = (a .+ b) ./ 2
    Γf = QuadraticTriangle(P1, P2, P3, mid(P1, P2), mid(P2, P3), mid(P3, P1))
    F(y) = 1 + y[1] * y[2] - y[3]^2
    for x0 in ((0.3, 0.3, 0.05), (0.5, -0.2, 0.01))
        rc = rule(WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γf; x0)); degree = 6)
        rf = rule(WeightedDomain(SurfaceTriangle(P1, P2, P3), InverseDistance(x0)); degree = 6)
        @test integrate(ξ -> F(Γf(ξ)), rc) ≈ integrate(F, rf) rtol = 1e-14
    end
end
