using CubatureRules, Test
import CubatureRules: MonicRecurrence, MomentWeight, OrdinaryMoments, ModifiedChebyshev,
                      MomentBreakdownError, monic, monomial_recurrence, wheeler,
                      jacobi_recurrence, gauss_from_recurrence, with_bits, digits_to_bits,
                      candidates, PolynomialDegree, passed
const CR = CubatureRules

# --- the two weights used throughout ------------------------------------------------------
# Lebesgue measure on [0,1] by its ordinary moments: the classically ill-conditioned case.
const LEBESGUE = OrdinaryMoments((k, T) -> one(T) / (k + 1); label = "Lebesgue on [0,1]")
const LEB_DOM = WeightedDomain(Interval(0, 1), LEBESGUE)

# Monic shifted Legendre on [0,1], the natural auxiliary family there.
const SHIFTED_LEGENDRE = MonicRecurrence((k, T) -> one(T) / 2,
                                         (k, T) -> k == 0 ? zero(T) : T(k)^2 / T(4 * (4 * k^2 - 1)))
# w(x) = log(1/x) on [0,1], which has no classical Gauss family. Modified moments against the
# monic shifted Legendre family: 1 at k = 0, else (-1)^k / (k(k+1)) over the leading
# coefficient binomial(2k, k).
const LOGW = MomentWeight(SHIFTED_LEGENDRE,
                          (k, T) -> k == 0 ? one(T) :
                                    T((-1)^k) / (T(k) * (k + 1) * T(binomial(big(2k), big(k))));
                          label = "log(1/x) on [0,1]")
const LOG_DOM = WeightedDomain(Interval(0, 1), LOGW)

# The two rules the assertions below keep coming back to, at 40 points and 50 digits: one
# from ordinary moments and one from a well-chosen auxiliary family. Built once each.
const HARD = rule(WeightedDomain(Interval(0, 1), LEBESGUE); degree = 79, digits = 50)
const EASY = rule(WeightedDomain(Interval(0, 1), LOGW); degree = 79, digits = 50)

"Gauss–Legendre on [0,1] at 1024 bits: the truth the ordinary-moment rules must reproduce."
function shifted_legendre_rule(n)
    return with_bits(1024) do
        x, w, _ = gauss_from_recurrence(n, jacobi_recurrence(0, 0), 1024)
        (sort([(xi + 1) / 2 for xi in x]), w ./ 2)
    end
end

@testset "Wheeler's algorithm" begin
    # Feeding a family its own modified moments (m = [1, 0, 0, …]) must return that family's
    # own coefficients: the sharpest available check on the a, b ≠ 0 path.
    self = MomentWeight(SHIFTED_LEGENDRE, (k, T) -> k == 0 ? one(T) : zero(T))
    with_bits(256) do
        α, β = wheeler(self, 20, BigFloat)
        @test all(a -> a == 1 // 2, α)
        @test all(k -> β[k] == SHIFTED_LEGENDRE.b(k - 1, BigFloat), 2:20)
        @test β[1] == 1
    end
    # Ordinary moments of Lebesgue measure on [0,1] give the shifted Legendre coefficients.
    with_bits(512) do
        α, β = wheeler(LEBESGUE, 12, BigFloat)
        @test maximum(abs, α .- 1 // 2) < 1e-100
        @test maximum(abs(β[k] - SHIFTED_LEGENDRE.b(k - 1, BigFloat)) for k in 2:12) < 1e-100
    end
    # monic() agrees with the hand-written shifted family, up to the interval it lives on
    m = monic(jacobi_recurrence(0, 0))
    @test m.b(3, Float64) ≈ jacobi_recurrence(0, 0).b(3, Float64)^2
    @test all(k -> monomial_recurrence().a(k, Float64) == 0, 0:5)

    # Too few moments to be a measure at all: Wheeler must say so rather than return numbers.
    bogus = OrdinaryMoments((k, T) -> k == 2 ? -one(T) : one(T))
    @test_throws MomentBreakdownError with_bits(128) do
        wheeler(bogus, 4, BigFloat)
    end
    @test_throws ArgumentError wheeler(LEBESGUE, 0, BigFloat)
end

@testset "selection" begin
    @test candidates(ModifiedChebyshev, LEB_DOM, PolynomialDegree(5)) == [ModifiedChebyshev()]
    @test family(rule(LEB_DOM; degree = 9)) == "ModifiedChebyshev"
    # a moment weight is stated on its own interval, so it is its own reference and the rule
    # is never built on [-1,1] and mapped
    @test CR.isreference(LEB_DOM)
    @test CR.reference(LEB_DOM) === LEB_DOM
    @test domain(rule(LEB_DOM; degree = 5)) === LEB_DOM
    @test measure(LEB_DOM) == 1
    @test npoints(rule(LEB_DOM; degree = 9)) == 5
    @test degree(rule(LEB_DOM; degree = 8)) == 9        # 5 points reach degree 9
    # exact arithmetic: the nodes are irrational, and the selector says so the same way it
    # does for every other Gauss family rather than failing inside the build
    @test_throws NoRuleError rule(LEB_DOM; degree = 5, T = Rational{BigInt})
end

@testset "ordinary moments: the ill-conditioned case" begin
    # The headline. In Float64 this construction collapses before n = 16 — β₁₃ goes negative
    # and the output stops describing a measure. Here the driver finds the precision the
    # conditioning demands, and the nodes come back right to the digits that were asked for.
    for (deg, n, digits) in ((9, 5, 30), (79, 40, 50))
        r = rule(LEB_DOM; degree = deg, digits = digits)
        @test npoints(r) == n
        xref, wref = shifted_legendre_rule(n)
        tol = exp10(-digits + 2)
        @test maximum(abs(BigFloat(nodes(r)[i][1]) - xref[i]) for i in 1:n) < tol
        @test maximum(abs(BigFloat(weights(r)[i]) - wref[i]) for i in 1:n) < tol
        @test all(>(0), weights(r))
        @test passed(check(r))
    end
    # The cost of the conditioning is reported, not hidden: 40 points from ordinary moments
    # loses well over fifty digits, and the certificate says so.
    c = certificate(HARD)
    @test c.cond > 1e40
    @test c.guard_digits > 100
    @test occursin("Wheeler", join(provenance(HARD).path, " "))
end

@testset "modified moments: a weight with no classical family" begin
    # w(x) = log(1/x) on [0,1]. Checked against ∫₀¹ xʲ log(1/x) dx = 1/(j+1)², which is exact
    # and is not what the rule was built from.
    for (n, digits) in ((5, 30), (40, 50))
        r = rule(LOG_DOM; degree = 2n - 1, digits = digits)
        @test npoints(r) == n
        @test all(>(0), weights(r))
        @test all(x -> 0 < x[1] < 1, nodes(r))
        worst = with_bits(digits_to_bits(digits) + 128) do
            maximum(0:(2n - 1)) do j
                got = sum(BigFloat(weights(r)[i]) * BigFloat(nodes(r)[i][1])^j for i in 1:n)
                abs(got - one(BigFloat) / BigFloat(j + 1)^2) * BigFloat(j + 1)^2
            end
        end
        @test worst < exp10(-digits + 2)
        @test passed(check(r))
    end
    # A well-chosen auxiliary family costs almost nothing: the same n that needs 150+ extra
    # digits from ordinary moments needs none here. This contrast is the point of the file.
    @test certificate(EASY).cond < 1e3
    @test certificate(HARD).cond > 1e40
    @test certificate(EASY).guard_digits < certificate(HARD).guard_digits ÷ 4
end

@testset "moments that cannot be met" begin
    # Moments of a signed measure: no positive rule exists, so no amount of precision helps
    # and the failure must be reported rather than papered over.
    signed = OrdinaryMoments((k, T) -> iseven(k) ? one(T) / (k + 1) : -one(T))
    @test_throws CR.RefinementError rule(WeightedDomain(Interval(0, 1), signed); degree = 11)
end
