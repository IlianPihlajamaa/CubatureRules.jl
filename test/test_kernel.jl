using CubatureRules, Test
const CR = CubatureRules

# Cauchy principal values and Hadamard finite parts with a Jacobi weight. The Hilbert
# transform of the weight is checked against classical closed forms and against tanh-sinh on
# the subtracted integrand; the rules are checked by `verify` and against the same reference.

# ⨍ or ⨎ of eˣ (b − x)^α (x − a)^β / (x − t)^order by subtraction: g = eˣ w is smooth at t,
# and what is left after removing its Taylor polynomial is integrated by tanh-sinh, with the
# limit value at a node that falls on t.
function kernel_reference(order, α, β, t, a, b; digits = 90)
    setprecision(BigFloat, 4 * digits) do
        A, B, T = BigFloat(a), BigFloat(b), BigFloat(t)
        α, β = BigFloat(α), BigFloat(β)
        w(x) = (B - x)^α * (x - A)^β
        L(x) = -α / (B - x) + β / (x - A)
        dL(x) = -α / (B - x)^2 - β / (x - A)^2
        g(x) = exp(x) * w(x)
        gt = g(T)
        dgt = gt * (1 + L(T))
        d2gt = gt * (1 + 2L(T) + L(T)^2 + dL(T))
        F = order == 1 ? (x -> x == T ? dgt : (g(x) - gt) / (x - T)) :
                         (x -> x == T ? d2gt / 2 : (g(x) - gt - dgt * (x - T)) / (x - T)^2)
        smooth = integrate(F, rule(TanhSinh(10), Interval(A, B); digits = digits + 10))
        order == 1 && return smooth + gt * log((B - T) / (T - A))
        return smooth + gt * (-1 / (B - T) - 1 / (T - A)) + dgt * log((B - T) / (T - A))
    end
end

@testset "Hilbert transform of a Jacobi weight" begin
    setprecision(BigFloat, 256) do
        for t in big.((-0.9, -0.4, 0.0, 0.25, 0.7))
            # Legendre: log((1 − t)/(1 + t)) and its derivative −2/(1 − t²)
            ρ, σ = CR.jacobi_hilbert(0, 0, t)
            @test abs(ρ - log((1 - t) / (1 + t))) < 1e-70
            @test abs(σ + 2 / (1 - t^2)) < 1e-70
            # Chebyshev: ⨍ (1 − x²)^{−1/2}/(x − t) = 0 and ⨍ (1 − x²)^{1/2}/(x − t) = −πt
            ρ, σ = CR.jacobi_hilbert(-1 // 2, -1 // 2, t)
            @test abs(ρ) < 1e-70 && abs(σ) < 1e-70
            ρ, σ = CR.jacobi_hilbert(1 // 2, 1 // 2, t)
            @test abs(ρ + big(π) * t) < 1e-70 && abs(σ + big(π)) < 1e-70
        end
    end
    # every branch of the series formulas, through the rules' integrals of eˣ against the
    # reference: non-integer, integer and nearly integer exponents, with t near either end
    # and in the middle
    for (α, β) in ((1 // 3, -1 // 2), (2, 0), (0, 3), (1, 1), (2, -1 // 2), (-1 // 2, 1), (1 // 10^9, 0)),
        t in (-8 // 10, -4 // 10, 1 // 10, 1 // 2, 9 // 10)
        dom = WeightedDomain(Interval(-1, 1), FinitePart(t; α, β))
        r = rule(dom; degree = 24, digits = 30)
        @test passed(verify(r))
        @test abs(integrate(exp, r) - kernel_reference(2, α, β, t, -1, 1)) < 1e-24
    end
end

@testset "principal-value and finite-part rules" begin
    cases = [(3 // 10, 0, 0, -1, 1), (9 // 10, 1 // 2, -1 // 2, -1, 1), (-99 // 100, -1 // 3, 2, -1, 1),
             (1 // 7, 1, 3 // 4, 0, 2), (1 // 3, -1 // 2, -1 // 2, -1, 1), (5 // 2, 0, 1 // 4, 2, 5)]
    for order in (1, 2), (t, α, β, a, b) in cases
        k = order == 1 ? PrincipalValue(t; α, β) : FinitePart(t; α, β)
        dom = WeightedDomain(Interval(a, b), k)
        r = rule(dom; degree = 40, digits = 40)
        @test family(r) == "SingularGauss"
        @test degree(r) >= 40 && npoints(r) <= degree(r) + 1
        @test any(x -> abs(x - t) < 1e-39, nodes(r))        # the singular point is a node
        @test all(x -> a < x < b, nodes(r))
        v = verify(r)
        @test passed(v) && v.sharp === true
        ref = kernel_reference(order, α, β, t, a, b)
        @test abs(integrate(exp, r) - ref) < 1e-37 * max(1, abs(ref))
    end
    # in Float64 and Float32, verified at the rules' own precision
    for T in (Float64, Float32), k in (PrincipalValue(0.3), FinitePart(0.3; α = 1 // 2))
        r = rule(WeightedDomain(Interval(-1, 1), k); degree = 12, T)
        @test eltype(r) == T && passed(verify(r))
        @test T(0.3) in nodes(r)
    end
    @test measure(WeightedDomain(Interval(-1, 1), PrincipalValue(1 // 3))) ≈ -log(2)
    @test measure(WeightedDomain(Interval(-1, 1), FinitePart(1 // 3))) ≈ -9 / 4
end

@testset "how the rules are built" begin
    dom(k) = WeightedDomain(Interval(-1, 1), k)
    # Hunter's rule: n Gauss nodes and t, degree 2n
    r = rule(dom(PrincipalValue(3 // 10)); degree = 20, digits = 30)
    @test npoints(r) == 11 && degree(r) == 20
    @test occursin("Gauss–Jacobi nodes and t", provenance(r).path[4])
    @test any(c -> c.key == "Hunter1972", provenance(r).citations)
    # t on a Gauss node of the first choice: one more point instead
    r = rule(dom(PrincipalValue(0)); degree = 30, digits = 30)
    @test npoints(r) == 17 && degree(r) == 32 && passed(verify(r))
    # the finite part with n + 1 nodes and degree 2n, centred: every odd polynomial is
    # integrated by symmetry, and the claim says so
    r = rule(dom(FinitePart(0)); degree = 30, digits = 30)
    @test degree(r) == 33 && passed(verify(r))
    @test occursin("zeros of ρ", provenance(r).path[4])
    # near the end, no n ≤ 2n₀ keeps every finite-part node inside: interpolatory instead
    r = rule(dom(FinitePart(999 // 1000)); degree = 44, digits = 30)
    @test npoints(r) == 45 && degree(r) == 44 && passed(verify(r))
    @test occursin("interpolatory", provenance(r).path[4])
    @test abs(integrate(exp, r) - kernel_reference(2, 0, 0, 999 // 1000, -1, 1)) < 1e-27
    # degree 0: the single node t
    r = rule(dom(FinitePart(3 // 10)); degree = 0, digits = 30)
    @test npoints(r) == 1 && abs(only(nodes(r)) - 3 // 10) < 1e-29 && passed(verify(r))
    # t within 1e-20 of an end
    r = rule(dom(FinitePart(1 - big(10)^-20)); degree = 20, digits = 30)
    @test passed(verify(r))
end

@testset "kernels: arguments and display" begin
    @test sprint(show, PrincipalValue(1 // 3)) == "PrincipalValue(1//3)"
    @test sprint(show, FinitePart(0.5; α = 1 // 2, β = 0)) == "FinitePart(0.5; α = 1//2, β = 0)"
    @test_throws ArgumentError PrincipalValue(0.3; α = -1)
    @test_throws ArgumentError FinitePart(Inf)
    @test_throws ArgumentError rule(WeightedDomain(Interval(0, 1), PrincipalValue(2)); degree = 4)
    @test_throws ArgumentError rule(WeightedDomain(Interval(0, 1), FinitePart(1)); degree = 4)
    @test_throws NoRuleError rule(WeightedDomain(Interval(-1, 1), PrincipalValue(0)); degree = 4, T = Rational{BigInt})
    @test [c.family for c in available(WeightedDomain(Interval(-1, 1), PrincipalValue(0.3)))] == ["SingularGauss"]
end
