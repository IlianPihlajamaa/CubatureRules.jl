using CubatureRules, LinearAlgebra, Test
const CR = CubatureRules

# The weakly singular volume kernel 1/|y − x₀| on a tetrahedron (DuffyCone). Besides `check`,
# which integrates the tetrahedral Dubiner polynomials over the cones from x₀ with its own
# polar quadrature on each face, an integral in spherical coordinates about x₀, which shares
# nothing with the cones: ∫_K 1/r dy = ∫_{S²} R(ω)²/2 dω, R(ω) the distance to the boundary.

K = Simplex((0, 0, 0), (1, 0, 0), (3 // 10, 9 // 10, 0), (1 // 5, 1 // 4, 4 // 5))

@testset "rules: exact, sharp, positive, interior" begin
    # `check` integrates the tetrahedral Dubiner polynomials over every face with a polar
    # quadrature refined to 128 bits: 20 s to a minute a rule, so it runs on two of them, and
    # the rest are held to the weight sum, ∫_K 1/r dy over the cones
    for (x0, d) in (((3 // 10, 3 // 10, 1 // 5), 3), ((0, 0, 0), 5))
        v = check(rule(WeightedDomain(K, InverseDistance(x0)); degree = d))
        @test passed(v) && v.exact && v.sharp === true && v.positive && v.interior
    end
    for x0 in ((3 // 10, 3 // 10, 1 // 5), (2 // 5, 3 // 10, 0), (1 // 2, 0, 0), (0, 0, 0), (3 // 10, 3 // 10, 1 // 1000))
        dom = WeightedDomain(K, InverseDistance(x0))
        μ = Float64(CR.measure(dom))
        for d in (0, 3, 5)
            r = rule(dom; degree = d)
            @test family(r) == "DuffyCone" && degree(r) >= d
            @test all(>(0), weights(r)) && all(y -> isinterior(y, K), nodes(r))
            @test sum(weights(r)) ≈ μ rtol = 1e-14
        end
    end
    r = rule(WeightedDomain(K, InverseDistance((3 // 10, 3 // 10, 1 // 5))); degree = 4, digits = 30)
    @test eltype(weights(r)) == BigFloat
    @test abs(sum(weights(r)) - setprecision(() -> CR.measure(r.domain), BigFloat, 256)) < 1e-28
    # x₀ outside: no rule (yet)
    @test_throws NoRuleError rule(WeightedDomain(K, InverseDistance((2, 2, 2))); degree = 3)
end

@testset "against spherical coordinates about x₀" begin
    x0 = [3 / 10, 3 / 10, 1 / 5]
    V = [Float64.(collect(v)) for v in CR.vertices(K)]
    # each face as (outward unit normal, a point on it)
    faces = map(((2, 3, 4), (1, 3, 4), (1, 2, 4), (1, 2, 3))) do (i, j, k)
        n = normalize(cross(V[j] - V[i], V[k] - V[i]))
        o = setdiff(1:4, (i, j, k))[1]
        dot(n, V[o] - V[i]) > 0 && (n = -n)
        (n, V[i])
    end
    R(ω) = minimum(dot(n, p - x0) / dot(n, ω) for (n, p) in faces if dot(n, ω) > 0)
    # R(ω) has kinks where the face changes, so the sphere rule converges slowly and erratically:
    # it recovers the volume, ∫ R³/3 dω, to about 1e-5, and ∫ R²/2 dω oscillates by as much
    # around the value the cones give (to 7e-7 at degree 1600)
    sphere = rule(Sphere{3}(); degree = 400)
    vol = abs(det(hcat(V[2] - V[1], V[3] - V[1], V[4] - V[1]))) / 6
    @test isapprox(sum(w * R(ω)^3 / 3 for (ω, w) in zip(nodes(sphere), weights(sphere))), vol; rtol = 1e-4)
    ref = sum(w * R(ω)^2 / 2 for (ω, w) in zip(nodes(sphere), weights(sphere)))
    r = rule(WeightedDomain(K, InverseDistance((3 // 10, 3 // 10, 1 // 5))); degree = 3)
    @test isapprox(sum(weights(r)), ref; rtol = 1e-4)
end
