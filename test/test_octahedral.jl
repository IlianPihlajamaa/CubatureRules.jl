using CubatureRules, Test, LinearAlgebra, Random
const CR = CubatureRules

# the 48 signed permutations, written out independently of the implementation
const OCT_PERMS = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
oct_images(x) = [[s1 * x[p[1]], s2 * x[p[2]], s3 * x[p[3]]]
                 for p in OCT_PERMS for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
p4of(x) = x[1]^4 + x[2]^4 + x[3]^4
p6of(x) = (x[1] * x[2] * x[3])^2

const ORBIT_PARAMS = Dict(:a1 => Float64[], :a2 => Float64[], :a3 => Float64[],
                          :b => [0.3], :c => [0.4], :d => [0.3, 0.5])

@testset "octahedral orbits" begin
    for (k, n) in ((:a1, 6), (:a2, 12), (:a3, 8), (:b, 24), (:c, 24), (:d, 48))
        o = CR.OctahedralOrbit(k)
        @test CR.orbit_size(o) == n
        pts = CR.orbit_points(o, ORBIT_PARAMS[k], Float64)
        @test length(pts) == n                       # no coincidences, no duplicates
        @test length(unique(pts)) == n
        @test all(x -> abs(sum(abs2, x) - 1) < 1e-14, pts)
        # the point set is invariant under the group, which is what makes it an orbit
        key(v) = round.(Float64.(v) .+ 0.0; digits = 12)    # `+ 0.0` folds -0.0 onto 0.0
        seen = Set(key(x) for x in pts)
        @test all(y -> key(y) in seen, Iterators.flatten(oct_images(x) for x in pts))
    end
    @test CR.nparams(CR.OctahedralOrbit(:a1)) == 0
    @test CR.nparams(CR.OctahedralOrbit(:b)) == 1
    @test CR.nparams(CR.OctahedralOrbit(:d)) == 2
    @test CR.nunknowns(CR.OctahedralOrbit(:d)) == 3
    @test_throws ArgumentError CR.OctahedralOrbit(:e)
    @test sprint(show, CR.OctahedralOrbit(:b)) == "OctahedralOrbit(:b)"
    # a zero coordinate is never negated: that would double the orbit
    @test length(CR.signed_orbit(CR.SVector(1.0, 0.0, 0.0))) == 6
    # parameters that leave the sphere are rejected rather than silently complexified
    @test_throws DomainError CR.orbit_representative(CR.OctahedralOrbit(:b), [0.8], Float64)
    @test_throws DomainError CR.orbit_representative(CR.OctahedralOrbit(:d), [0.9, 0.9], Float64)
end

@testset "octahedral structures" begin
    s = CR.OctahedralStructure([:a1, :a3, :b, :d])
    @test CR.npoints(s) == 6 + 8 + 24 + 48
    @test CR.nunknowns(s) == 1 + 1 + 2 + 3
    @test CR.param_offsets(s) == [0, 1, 2, 4]
    @test occursin("86 points", sprint(show, s))
    xs, ws = CR.expand(s, [0.01, 0.02, 0.03, 0.3, 0.04, 0.3, 0.5])
    @test length(xs) == length(ws) == 86
    @test all(x -> abs(sum(abs2, x) - 1) < 1e-14, xs)
    @test sort(unique(ws)) == [0.01, 0.02, 0.03, 0.04]
    @test_throws ArgumentError CR.expand(s, [0.01, 0.02])
end

@testset "octahedral invariants" begin
    # an invariant is constant on an orbit — the fact that makes the system cheap
    for (k, par) in ORBIT_PARAMS
        o = CR.OctahedralOrbit(k)
        pts = CR.orbit_points(o, par, Float64)
        P4, P6, _, _ = CR.orbit_invariants(o, par, Float64)
        @test maximum(x -> abs(p4of(x) - P4) + abs(p6of(x) - P6), pts) < 1e-15
    end
    # the basis {p₄^a p₆^b : 4a+6b ≤ n} is independent and spans the invariants. Checked
    # against the group average of every monomial of degree ≤ n, not assumed.
    Random.seed!(20260922)
    for n in (4, 6, 8, 12, 16)
        exps = CR.invariant_exponents(n)
        pts = [(v = randn(3); v ./ norm(v)) for _ in 1:400]
        A = [p4of(x)^a * p6of(x)^b for x in pts, (a, b) in exps]
        mons = [(i, j, k) for i in 0:n for j in 0:n for k in 0:n if i + j + k <= n]
        B = [sum(y -> y[1]^α[1] * y[2]^α[2] * y[3]^α[3], oct_images(x)) / 48 for x in pts, α in mons]
        rk(M) = (s = svdvals(M); count(>(1e-9 * s[1]), s))
        @test rk(A) == length(exps)          # independent
        @test rk(B) == length(exps)          # and spanning: nothing else is invariant
    end
    @test CR.invariant_exponents(3) == [(0, 0)]
    @test CR.invariant_exponents(6) == [(0, 0), (1, 0), (0, 1)]
    # the right-hand sides, against a product rule of ample degree
    for (a, b) in CR.invariant_exponents(12)
        r = rule(Sphere{3}(); degree = 4a + 6b + 2)
        num = integrate(x -> p4of(x)^a * p6of(x)^b, r)
        @test Float64(CR.invariant_moment(a, b)) ≈ num rtol = 1e-12
    end
end

@testset "octahedral moment system" begin
    s = CR.OctahedralStructure([:a1, :b, :d])
    sys = CR.OctahedralMomentSystem(s, 11, Float64)
    @test CR.n_equations(sys) == length(CR.invariant_exponents(11))
    @test CR.n_unknowns(sys) == 6
    θ = [0.05, 0.03, 0.35, 0.04, 0.3, 0.5]
    r, J = sys(θ)
    @test size(J) == (CR.n_equations(sys), 6)
    # the analytic Jacobian against a central difference
    for j in eachindex(θ)
        h = 1e-7
        θp = copy(θ); θp[j] += h
        θm = copy(θ); θm[j] -= h
        @test (sys(θp)[1] - sys(θm)[1]) / 2h ≈ J[:, j] rtol = 1e-5
    end
    # the residual is the rule's error on each invariant, so a rule that integrates them
    # exactly has residual zero — the degree-5 case can be solved by hand
    s5 = CR.OctahedralStructure([:a1, :a3])
    sys5 = CR.OctahedralMomentSystem(s5, 5, Float64)
    w = sys5([0.0, 0.0])[2] \ Float64[CR.invariant_moment(0, 0), CR.invariant_moment(1, 0)]
    @test w[1] / (4π) ≈ 1 / 15 rtol = 1e-13          # the classical Lebedev 14-point weights,
    @test w[2] / (4π) ≈ 3 / 40 rtol = 1e-13          # derived here from the equations alone
    @test maximum(abs, sys5(w)[1]) < 1e-14
end

@testset "octahedral refinement" begin
    s5 = CR.OctahedralStructure([:a1, :a3])
    sys5 = CR.OctahedralMomentSystem(s5, 5, Float64)
    w = sys5([0.0, 0.0])[2] \ Float64[CR.invariant_moment(0, 0), CR.invariant_moment(1, 0)]
    for digits in (30, 80)
        bits = CR.digits_to_bits(digits)
        θ, res, guard = CR.refine_octahedral(s5, 5, Float64[w[1], w[2]], bits)
        @test res.converged
        @test guard > 0
        xs, ws = CR.with_bits(() -> CR.expand(s5, θ), bits + guard)
        r5 = CR.QuadratureRule(xs, ws, Sphere{3}(), PolynomialDegree(5),
                               CR.Provenance(family = "test", derivation = CR.Derived()))
        v = check(r5)
        @test v.exact && v.sharp === true && v.positive && v.interior && v.weights_sum_ok
        # the sum has to be formed at the rule's own precision: BigFloat arithmetic uses the
        # ambient setting, and the default 256 bits is coarser than an 80-digit rule
        dev = CR.with_bits(() -> abs(sum(ws) - 4 * BigFloat(π)), bits + guard)
        @test dev < big(10.0)^(-digits + 2)
    end
end
