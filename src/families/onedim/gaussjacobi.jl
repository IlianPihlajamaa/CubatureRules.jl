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

"Guard bits for Newton on a three-term recurrence: O(log n), which is right there (§2.7)."
gj_guard_bits(n) = 24 + 2 * ceil(Int, log2(n + 1))

"""
    gauss_jacobi_work(n, α, β, bits) -> (x, w, iterations)

`n`-point Gauss–Jacobi nodes (ascending) and weights as BigFloats of precision `bits`.
Callers choose `bits` including any guard; nothing here reads the ambient precision.
"""
function gauss_jacobi_work(n::Integer, α, β, bits::Integer)
    n >= 1 || throw(ArgumentError("need at least one node"))
    # --- seed: Golub–Welsch in Float64
    a64 = [Float64(jacobi_a(k, float(α), float(β))) for k in 0:(n - 1)]
    b64 = [Float64(jacobi_b(k, float(α), float(β))) for k in 1:(n - 1)]
    seed = n == 1 ? a64 : eigvals(SymTridiagonal(a64, b64))
    sort!(seed)
    return with_bits(bits) do
        A, B = big(α), big(β)
        a = [jacobi_a(k, A, B) for k in 0:n]
        b = [jacobi_b(k, A, B) for k in 1:n]
        μ = BigFloat(jacobi_mass(α, β))
        x = BigFloat.(seed)
        tol = ldexp(BigFloat(1), -(bits - 6))
        iters = 0
        for i in 1:n
            xi = x[i]
            last_small = false
            for it in 1:100
                q, dq, _ = jacobi_eval(n, xi, a, b)
                δ = q / dq
                xi -= δ
                iters = max(iters, it)
                small = abs(δ) <= tol * max(abs(xi), one(xi))
                small && last_small && break
                last_small = small
                it == 100 && throw(RefinementError("GaussJacobi", "Newton on P_$n did not converge at node $i"))
            end
            x[i] = xi
        end
        w = similar(x)
        for i in 1:n
            _, _, s = jacobi_eval(n, x[i], a, b)
            w[i] = μ / s
        end
        x, w, iters
    end
end

# Unnormalised (q₀ = 1) orthonormal recurrence: returns q_n(x), q_n'(x) and Σ_{k<n} q_k(x)².
function jacobi_eval(n, x, a, b)
    q0, q1 = zero(x), one(x)
    d0, d1 = zero(x), zero(x)
    s = one(x)
    for k in 0:(n - 1)
        bk = k == 0 ? zero(x) : b[k]
        q2 = ((x - a[k + 1]) * q1 - bk * q0) / b[k + 1]
        d2 = (q1 + (x - a[k + 1]) * d1 - bk * d0) / b[k + 1]
        q0, q1 = q1, q2
        d0, d1 = d1, d2
        k < n - 1 && (s += q1^2)
    end
    return q1, d1, s
end

"Newton-correction residual `max |q_n(x̂)/q_n'(x̂)|` of rounded nodes, evaluated at `bits`."
function gj_residual(xhat, α, β, bits)
    n = length(xhat)
    with_bits(bits) do
        A, B = big(α), big(β)
        a = [jacobi_a(k, A, B) for k in 0:n]
        b = [jacobi_b(k, A, B) for k in 1:n]
        maximum(xhat) do xi
            q, dq, _ = jacobi_eval(n, BigFloat(xi), a, b)
            abs(q / dq)
        end
    end
end

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
