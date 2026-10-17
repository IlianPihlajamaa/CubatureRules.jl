using CubatureRules, LinearAlgebra, StaticArrays, Test
const CR = CubatureRules

# KernelRuleBuilder: the boundary-element kernel rules built fast in Float64 (src/fast). Against
# the certified rules of DuffyGauss, DuffySinh and DuffyFinitePart, and through `check`; the
# fallback outside the validated range; no allocations once warm.

T = SurfaceTriangle((0.1, 0.2, 0.3), (1.3, -0.4, 0.7), (0.2, 1.1, -0.5))
V = CR.vertices(T)
n = normalize(cross(V[2] - V[1], V[3] - V[1]))
c = (V[1] + V[2] + V[3]) / 3
points = [("vertex", V[2]), ("edge", (V[1] + V[3]) / 2), ("inside", c),
          ("above", c + 1e-3 * n), ("below, close", c - 1e-9 * n), ("beside", (V[1] + V[2]) / 2 + 0.05 * normalize(V[1] + V[2] - 2V[3])),
          ("in the plane, outside", V[1] + 0.1 * normalize(V[1] - c))]
kernel(K, x0) = K === InverseDistanceGradient ? K(Tuple(x0), (0.6, 0.0, 0.8)) : K(Tuple(x0))

# the integrals of a few polynomials, and how much they may differ for a Float64 rule
fs = (y -> 1.0, y -> 1 + y[1] - 2y[2] * y[3], y -> y[1]^3 - y[2]^2 * y[3] + 0.5)

@testset "against the certified rules" begin
    for K in (InverseDistance, InverseDistanceCubed, InverseDistanceGradient), d in (3, 8)
        b = KernelRuleBuilder(K, d)
        for (name, x0) in points
            k = kernel(K, x0)
            x, _ = b(T, k)
            fast = [integrate(f, b, T, k) for f in fs]
            r = rule(WeightedDomain(T, k); degree = d)
            ref = [integrate(f, r) for f in fs]
            scale = maximum(f -> sum(abs.(weights(r)) .* abs.(f.(nodes(r)))), fs)
            @test b.place !== :fallback
            @test maximum(abs.(fast .- ref)) <= 1e-13 * scale
            b.place !== :near && @test length(x) == npoints(r)
        end
    end
end

@testset "check" begin
    for (K, d, x0) in ((InverseDistance, 5, c), (InverseDistance, 5, c + 1e-3 * n),
                       (InverseDistanceCubed, 4, (V[1] + V[3]) / 2), (InverseDistanceGradient, 4, c + 1e-2 * n))
        r = rule(KernelRuleBuilder(K, d), T, kernel(K, x0))
        v = check(r)
        @test passed(v) && v.interior
        @test occursin("(fast)", family(r))
    end
end

@testset "x₀ very close to a triangle far from the origin" begin
    # the height of x₀ above the plane is computed in Float64x2: in Float64, n·(X − V₁) would be
    # off by eps |X|, a relative error of 1e-4 here, which ∫ 1/r³ ≈ Ω/h inherits
    S = SurfaceTriangle((3.1, 2.2, 1.3), (3.4, 2.15, 1.35), (3.2, 2.45, 1.25))
    W = CR.vertices(S)
    m = normalize(cross(W[2] - W[1], W[3] - W[1]))
    x0 = (W[1] + W[2] + W[3]) / 3 + 1e-12 * m
    b = KernelRuleBuilder(InverseDistanceCubed, 3)
    k = InverseDistanceCubed(Tuple(x0))
    ref = integrate(y -> 1.0, rule(WeightedDomain(S, k); degree = 3, digits = 20))
    @test integrate(y -> 1.0, b, S, k) ≈ ref rtol = 1e-11
end

@testset "x₀ just above a thin triangle, close to two long edges" begin
    # the far edges' singularities are η/κ² ≈ 170 times the integral; without that in the
    # angular line count, the gradient here was off by 4.5e-12
    S = SurfaceTriangle((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.45, 0.012, 0.0))
    k = InverseDistanceGradient((0.45, 0.006, 1e-9), Tuple(normalize([1.0, 0.3, 0.01])))
    b = KernelRuleBuilder(InverseDistanceGradient, 0)
    ref = integrate(y -> 1.0, rule(WeightedDomain(S, k); degree = 0, digits = 20))
    @test abs(integrate(y -> 1.0, b, S, k) - ref) < 1e-13
    @test b.place === :near
end

@testset "triangles in the plane" begin
    P = Simplex((0.0, 0.0), (1.0, 0.0), (0.3, 0.8))
    for K in (InverseDistance, InverseDistanceCubed), x0 in ((0.4, 0.3), (0.5, -0.02), (0.0, 0.0))
        b = KernelRuleBuilder(K, 6; dim = 2)
        k = K((x0))
        ref = rule(WeightedDomain(P, k); degree = 6)
        @test integrate(y -> 1 + y[1] * y[2], b, P, k) ≈ integrate(y -> 1 + y[1] * y[2], ref) rtol = 1e-12
    end
end

@testset "outside the validated range: the certified construction" begin
    # x₀ on the triangle within 1e-11 of an edge: a piece thinner than the sweep covered
    b = KernelRuleBuilder(InverseDistance, 5)
    x0 = (V[1] + V[2]) / 2 + 1e-11 * normalize(c - (V[1] + V[2]) / 2)
    k = InverseDistance(Tuple(x0))
    x, w = b(T, k)
    @test b.place === :fallback && b.fallbacks == 1
    r = rule(WeightedDomain(T, k); degree = 5)
    @test x == nodes(r) && w == weights(r)
end

# measured inside a function: at a testset's scope the returned pair of vectors escapes and is
# boxed (32 bytes on Julia 1.11), which code that uses the rule does not pay
first_coordinate(y) = y[1]
function allocations(b, T, k)
    b(T, k)
    integrate(first_coordinate, b, T, k)
    return @allocated(b(T, k)) + @allocated(integrate(first_coordinate, b, T, k))
end

@testset "buffers, allocations, arguments" begin
    b = KernelRuleBuilder(InverseDistanceGradient, 7)
    k1, k2 = kernel(InverseDistanceGradient, c), kernel(InverseDistanceGradient, c + 1e-4 * n)
    @test allocations(b, T, k1) == 0
    @test allocations(b, T, k2) == 0
    @test allocations(KernelRuleBuilder(InverseDistance, 9), T, InverseDistance(Tuple(c))) == 0
    @test_throws ArgumentError b(T, InverseDistance(Tuple(c)))
    @test_throws ArgumentError KernelRuleBuilder(InverseDistance, 20)
    @test_throws ArgumentError KernelRuleBuilder(Float64, 3)
    v = CR.validation(b)
    @test v === nothing || (v.geometries > 0 && v.worst < 1e-8)
end
