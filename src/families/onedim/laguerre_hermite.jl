# Gauss–Laguerre and Gauss–Hermite (PLAN §6 Tier 1): the classical rules on unbounded
# domains. Neither is in the upstream packages we delegate to, and both are three lines of
# recurrence coefficients on top of the shared Gauss driver.
#
#   Laguerre, weight x^α e^{-x} on [0, ∞):   a_k = 2k + α + 1,  b_k = √(k (k+α)),  μ₀ = Γ(α+1)
#   Hermite,  weight e^{-x²} on ℝ:           a_k = 0,           b_k = √(k/2),      μ₀ = √π

"""
    GaussLaguerre(α = 0)

Gauss–Laguerre rules for the weight `x^α e^{-x}` on `[0, ∞)`. An `n`-point rule has
polynomial degree `2n - 1` against that weight, and lives on [`LaguerreRay`](@ref)`(α)`.
"""
struct GaussLaguerre{T} <: RuleFamily
    α::T
    function GaussLaguerre(α::T = 0) where {T}
        α > -1 || throw(ArgumentError("Gauss–Laguerre needs α > -1, got $α"))
        return new{T}(α)
    end
end

"""
    GaussHermite()

Gauss–Hermite rules for the weight `e^{-x²}` on the whole line. An `n`-point rule has
polynomial degree `2n - 1` against that weight, and lives on [`HermiteLine`](@ref)`()`.
"""
struct GaussHermite <: RuleFamily end

derivation(::Type{<:GaussLaguerre}) = Derived()
derivation(::Type{GaussHermite}) = Derived()
describe_family(f::GaussLaguerre) = iszero(f.α) ? "GaussLaguerre" : "GaussLaguerre($(f.α))"
family_name(::GaussLaguerre) = "GaussLaguerre"

const SZEGO_LAGUERRE = SZEGO_1939

"""
    laguerre_recurrence(α)
    hermite_recurrence()

The orthonormal recurrences of the two weights, for the shared Gauss driver.
"""
laguerre_recurrence(α) = Recurrence((k, T) -> T(2k) + T(α) + one(T), (k, T) -> sqrt(T(k) * (T(k) + T(α))),
                                    T -> T(laguerre_mass(α)))
hermite_recurrence() = Recurrence((k, T) -> zero(T), (k, T) -> sqrt(T(k) / 2), T -> sqrt(T(π)))

recurrence_of(f::GaussLaguerre) = laguerre_recurrence(f.α)
recurrence_of(::GaussHermite) = hermite_recurrence()
home_domain(f::GaussLaguerre) = LaguerreRay(f.α)
home_domain(::GaussHermite) = HermiteLine()

candidates(::Type{<:GaussLaguerre}, dom::LaguerreDomain, ::PolynomialDegree) =
    isreference(dom) ? [GaussLaguerre(dom.weight.α)] : GaussLaguerre[]
candidates(::Type{GaussHermite}, dom::HermiteDomain, ::PolynomialDegree) =
    isreference(dom) ? [GaussHermite()] : GaussHermite[]

const ClassicalUnbounded = Union{GaussLaguerre,GaussHermite}

npoints(::ClassicalUnbounded, dom, degree::Integer) = gj_points(degree)
claimed_degree(::ClassicalUnbounded, dom, degree) = 2gj_points(degree) - 1
degree_range(::ClassicalUnbounded, dom::WeightedDomain) = 0:typemax(Int)
degree_range(::ClassicalUnbounded, dom) = 1:0
degree_for_npoints(::ClassicalUnbounded, dom, n::Integer) = 2n - 1
properties(::GaussLaguerre, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)
properties(::GaussHermite, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = false)

function build(f::ClassicalUnbounded, dom::WeightedDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("$(describe_family(f)) nodes are irrational; $(T) is not supported"))
    n = gj_points(degree)
    rec = recurrence_of(f)
    guard = gj_guard_bits(n) + 8n ÷ 10          # Laguerre nodes spread to ~4n, so they need room
    x, w, iters = gauss_from_recurrence(n, rec, ctx.bits + guard; name = family_name(f))
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    res = gauss_residual(xs, rec, 2ctx.bits + guard)
    cert = Certificate(equations = "p_$n(x_i) = 0 for the orthonormal polynomials of the weight " *
                                   "(Newton correction |p/p'|)",
                       residual = BigFloat(res; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = 1.0, iterations = iters)
    prov = Provenance(family = family_name(f), derivation = Derived(),
                      path = ["seed: Golub–Welsch eigenvalues in Float64",
                              "refine: Newton on the three-term recurrence at $(ctx.bits + guard) bits",
                              "weights: Christoffel function"],
                      seed_source = "Golub–Welsch (Float64)", citations = [SZEGO_LAGUERRE, GOLUB_WELSCH_1969],
                      symmetry = properties(f, dom, degree).symmetry)
    return QuadratureRule(xs, ws, home_domain(f), PolynomialDegree(2n - 1), prov, cert)
end
