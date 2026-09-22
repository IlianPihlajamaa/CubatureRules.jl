# Gauss–Jacobi rules (PLAN §6 Tier 1), included in v0.1 because the conical product needs
# them. Seed: Golub–Welsch in Float64. Refine: Newton on the three-term recurrence in
# BigFloat. Weights from the Christoffel function w_i = 1 / Σ_k p_k(x_i)².

"""
    GaussJacobi(α, β)
    GaussJacobi()

Gauss–Jacobi rules for the weight `(1-x)^α (1+x)^β` on `[-1, 1]`; `GaussJacobi()` (and
[`GaussLegendre`](@ref)`()`) is the unweighted case. An `n`-point rule has polynomial
degree `2n - 1`.
"""
struct GaussJacobi{T} <: RuleFamily
    α::T
    β::T
    function GaussJacobi{T}(α, β) where {T}
        (α > -1 && β > -1) || throw(ArgumentError("Gauss–Jacobi needs α, β > -1, got ($α, $β)"))
        return new{T}(α, β)
    end
end
GaussJacobi(α::T, β::T) where {T} = GaussJacobi{T}(α, β)
GaussJacobi(α, β) = GaussJacobi(promote(α, β)...)
GaussJacobi() = GaussJacobi(0, 0)

"""
    GaussLegendre()

Gauss–Legendre rules: `GaussJacobi(0, 0)`.
"""
GaussLegendre() = GaussJacobi(0, 0)

derivation(::Type{<:GaussJacobi}) = Derived()
conical_compatible(::Type{<:GaussJacobi}) = true
describe_family(f::GaussJacobi) = iszero(f.α) && iszero(f.β) ? "GaussLegendre" : "GaussJacobi($(f.α), $(f.β))"
family_name(::GaussJacobi) = "GaussJacobi"

const JacobiInterval = WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight}

candidates(::Type{<:GaussJacobi}, dom::Interval, ::PolynomialDegree) = [GaussJacobi()]
candidates(::Type{<:GaussJacobi}, dom::JacobiInterval, ::PolynomialDegree) =
    isreference(dom.base) ? [GaussJacobi(dom.weight.α, dom.weight.β)] : GaussJacobi[]

gj_points(degree) = max(1, cld(degree + 1, 2))
npoints(::GaussJacobi, dom, degree::Integer) = gj_points(degree)
claimed_degree(::GaussJacobi, dom, degree) = 2gj_points(degree) - 1
degree_range(::GaussJacobi, dom) = 0:typemax(Int)
properties(::GaussJacobi, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)
degree_for_npoints(::GaussJacobi, dom, n::Integer) = 2n - 1

gj_domain(f::GaussJacobi) = iszero(f.α) && iszero(f.β) ? Interval() : WeightedDomain(Interval(), JacobiWeight(f.α, f.β))

const SZEGO_1939 = Citation(key = "Szego1939", authors = ["Gábor Szegő"],
                            title = "Orthogonal Polynomials", journal = "American Mathematical Society Colloquium Publications",
                            year = 1939, volume = "23")
const GOLUB_WELSCH_1969 = Citation(key = "GolubWelsch1969", authors = ["Gene H. Golub", "John H. Welsch"],
                                   title = "Calculation of Gauss quadrature rules",
                                   journal = "Mathematics of Computation", year = 1969, volume = "23",
                                   pages = "221--230", doi = "10.1090/S0025-5718-69-99647-1")

"""
    jacobi_recurrence(α, β)

The orthonormal recurrence of the Jacobi weight `(1-x)^α (1+x)^β`, in the form the shared
Gauss driver takes.
"""
jacobi_recurrence(α, β) = Recurrence((k, T) -> jacobi_a(k, T(α), T(β)), (k, T) -> jacobi_b(k, T(α), T(β)),
                                     T -> T(jacobi_mass(α, β)))

"""
    gauss_jacobi_work(n, α, β, bits) -> (x, w, iterations)

`n`-point Gauss–Jacobi nodes (ascending) and weights as BigFloats of precision `bits`.
"""
gauss_jacobi_work(n::Integer, α, β, bits::Integer) =
    gauss_from_recurrence(n, jacobi_recurrence(α, β), bits; name = "GaussJacobi")

gj_residual(xhat, α, β, bits) = gauss_residual(xhat, jacobi_recurrence(α, β), bits)

function build(f::GaussJacobi, dom, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("Gauss–Jacobi nodes are irrational; $(T) is not supported"))
    n = gj_points(degree)
    guard = gj_guard_bits(n)
    x, w, iters = gauss_jacobi_work(n, f.α, f.β, ctx.bits + guard)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    res = gj_residual(xs, f.α, f.β, 2ctx.bits + guard)
    cert = Certificate(equations = "P_$n^($(f.α),$(f.β))(x_i) = 0 (Newton correction |P/P'|)",
                       residual = BigFloat(res; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = 1.0, iterations = iters)
    prov = Provenance(family = "GaussJacobi", derivation = Derived(),
                      path = ["seed: Golub–Welsch eigenvalues in Float64",
                              "refine: Newton on the three-term recurrence at $(ctx.bits + guard) bits",
                              "weights: Christoffel function"],
                      seed_source = "Golub–Welsch (Float64)",
                      citations = [GOLUB_WELSCH_1969, SZEGO_1939])
    return QuadratureRule(xs, ws, gj_domain(f), PolynomialDegree(2n - 1), prov, cert)
end
