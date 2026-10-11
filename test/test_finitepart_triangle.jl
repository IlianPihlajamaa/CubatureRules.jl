using CubatureRules, LinearAlgebra, Test
const CR = CubatureRules

# The finite parts, with respect to r = |y − x₀|, of the hypersingular kernel 1/r³ and the
# strongly singular (y − x₀)·e/r³ on a triangle with x₀ on it (DuffyFinitePart). Besides `check`,
# which compares with finite parts taken in polar coordinates about x₀, two closed forms and the
# additivity over triangles that makes the finite part the one boundary elements need.

T = Simplex((0, 0), (1, 0), (3 // 10, 8 // 10))
e = (3 // 5, 4 // 5)

"Foot distance and edge parameters of the edge (p, q) seen from x₀."
function edge_frame(x0, p, q)
    p, q, x0 = big.(collect(p)), big.(collect(q)), big.(collect(x0))
    u = (q - p) / norm(q - p)
    f = p + dot(x0 - p, u) * u
    return norm(x0 - f), dot(p - f, u), dot(q - f, u), u
end

@testset "rules: exact, sharp, interior, in the plane and in space" begin
    S3 = SurfaceTriangle((0, 0, 0), (1, 0, 1), (3 // 10, 8 // 10, 0))
    for (base, x0, ee) in ((T, (0, 0), e), (T, (6 // 10, 0), e), (T, (2 // 5, 3 // 10), e),
                           (S3, (13 // 30, 8 // 30, 1 // 3), (3 // 5, 0, 4 // 5)))
        for w in (InverseDistanceCubed(x0), InverseDistanceGradient(x0, ee)), d in (0, 3, 7)
            dom = WeightedDomain(base, w)
            r = rule(dom; degree = d)
            v = check(r)
            @test family(r) == "DuffyFinitePart" && npoints(r) == npoints(DuffyFinitePart(), dom, d)
            @test passed(v) && v.interior
        end
    end
    r = rule(WeightedDomain(T, InverseDistanceCubed((2 // 5, 3 // 10))); degree = 6, digits = 40)
    @test passed(check(r)) && eltype(weights(r)) == BigFloat
    # x₀ off the triangle: an ordinary (near-singular) integral, for DuffySinh
    @test family(rule(WeightedDomain(T, InverseDistanceCubed((2, 2))); degree = 3)) == "DuffySinh"
end

@testset "closed forms" begin
    V = CR.vertices(T)
    setprecision(BigFloat, 256) do
        # ⨎_T 1/r³ dy = −Σ (sin ψ_q − sin ψ_p)/h over the edges not through x₀, sin ψ = τ/√(h² + τ²)
        for x0 in ((0, 0), (6 // 10, 0), (2 // 5, 3 // 10))
            fp = -sum(((p, q) for (p, q) in ((V[1], V[2]), (V[2], V[3]), (V[3], V[1]))
                       if edge_frame(x0, p, q)[1] > 1e-30); init = big(0)) do (p, q)
                h, τp, τq, _ = edge_frame(x0, p, q)
                (τq / sqrt(h^2 + τq^2) - τp / sqrt(h^2 + τp^2)) / h
            end
            r = rule(WeightedDomain(T, InverseDistanceCubed(x0)); degree = 2, digits = 40)
            @test abs(sum(weights(r)) - fp) < 1e-36 * abs(fp)
        end
        # x₀ inside: PV ∫_T (y − x₀)·e/r³ dy = −∮ (n·e)/r ds (the divergence theorem; the circle
        # around x₀ contributes ∮ ê dθ = 0), with ∫ ds/r = asinh(τ_q/h) − asinh(τ_p/h) per edge
        x0 = (2 // 5, 3 // 10)
        pv = -sum(((V[1], V[2]), (V[2], V[3]), (V[3], V[1]))) do (p, q)
            h, τp, τq, u = edge_frame(x0, p, q)
            n = [u[2], -u[1]]                                 # outward: the vertices run counterclockwise
            dot(n, big.(collect(e))) * (asinh(τq / h) - asinh(τp / h))
        end
        r = rule(WeightedDomain(T, InverseDistanceGradient(x0, e)); degree = 2, digits = 40)
        @test abs(sum(weights(r)) - pv) < 1e-36 * abs(pv)
    end
end

@testset "additive over triangles" begin
    # T cut along the line from its first vertex through x₀, which then lies on the common edge
    x0 = (2 // 5, 3 // 10)
    P = (32 // 53, 24 // 53)                                  # where the line meets the opposite edge
    T1, T2 = Simplex((0, 0), (1, 0), P), Simplex((0, 0), P, (3 // 10, 8 // 10))
    f(y) = 1 + y[1] * y[2] - 2y[2]^3
    for w in (InverseDistanceCubed(x0), InverseDistanceGradient(x0, e))
        rT = rule(WeightedDomain(T, w); degree = 4, digits = 30)
        r1 = rule(WeightedDomain(T1, w); degree = 4, digits = 30)
        r2 = rule(WeightedDomain(T2, w); degree = 4, digits = 30)
        @test abs(integrate(f, rT) - integrate(f, r1) - integrate(f, r2)) < 1e-26
    end
end

@testset "near-singular: x₀ off the triangle" begin
    # above a flat triangle at the height H: ∫ 1/r³ = Ω/H and the double layer
    # ∫ (y − x₀)·n/r³ = −Ω, Ω the solid angle (Van Oosterom & Strackee), however close x₀ is
    V = ((0, 0, 0), (1, 0, 0), (3 // 10, 8 // 10, 0))
    S3 = SurfaceTriangle(V...)
    function solid_angle(x0)
        R = [big.(collect(v)) .- big.(collect(x0)) for v in V]
        r = norm.(R)
        num = abs(dot(R[1], cross(R[2], R[3])))
        den = r[1] * r[2] * r[3] + dot(R[1], R[2]) * r[3] + dot(R[1], R[3]) * r[2] + dot(R[2], R[3]) * r[1]
        return 2 * atan(num, den)
    end
    for foot in ((1 // 3, 1 // 4), (1 // 2, -1 // 10)), H in (1 // 10, 1 // 10^4, 1 // 10^9)
        x0 = (foot..., H)
        Ω = solid_angle(x0)
        r3 = rule(WeightedDomain(S3, InverseDistanceCubed(x0)); degree = 5)
        rg = rule(WeightedDomain(S3, InverseDistanceGradient(x0, (0, 0, 1))); degree = 5)
        @test family(r3) == "DuffySinh" && degree(r3) >= 5 && degree(rg) >= 5
        @test sum(weights(r3)) ≈ Ω / H rtol = 1e-14
        @test sum(weights(rg)) ≈ -Ω rtol = 1e-14
    end
    # against the polar reference, with e tilted, and in the plane outside the triangle
    for dom in (WeightedDomain(S3, InverseDistanceGradient((1 // 3, 1 // 4, 1 // 1000), (3 // 5, 0, 4 // 5))),
                WeightedDomain(T, InverseDistanceCubed((1 // 2, -1 // 100))))
        v = check(rule(dom; degree = 4))
        @test passed(v) && v.interior
    end
end
