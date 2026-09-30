# Gauss rules for a weight given as a function (PLAN §6 Tier 2: Gautschi's multiple-component
# discretisation). See domains/functionweight.jl for how the weight is discretised.
#
# The recurrence of the discretised measure comes from the discrete Stieltjes procedure,
# written here rather than taken from a package: PolyChaos.jl, the one Julia implementation,
# computes in Float64 whatever its input (β correct to 1e-16 on a 256-bit measure), and
# QuadGK's arbitrary-precision construction for a weight function evaluates Chebyshev series
# at every node on every step, 30× slower on these sizes. The procedure itself is a dozen
# lines. It can lose accuracy when the number of coefficients approaches the size of the
# discrete measure, so the discretisation keeps at least twice as many points as
# coefficients, and the agreement test between resolutions catches the rest.

"""
    discrete_stieltjes(N, x, λ) -> (α, β)

The first `N` monic recurrence coefficients of the discrete measure `Σᵢ λᵢ δ(x − xᵢ)`, by the
Stieltjes procedure: `β[1]` is the mass, and `π_{k+1} = (x − α_k) π_k − β_k π_{k−1}`.
"""
function discrete_stieltjes(N::Integer, x::AbstractVector{T}, λ::AbstractVector{T}) where {T}
    α, β = zeros(T, N), zeros(T, N)
    p0, p1 = zeros(T, length(x)), ones(T, length(x))
    s = sum(λ)
    β[1], α[1] = s, sum(λ .* x) / s
    for k in 1:(N - 1)
        p0, p1 = p1, (x .- α[k]) .* p1 .- (k == 1 ? zero(T) : β[k]) .* p0
        s1 = sum(λ .* p1 .^ 2)
        α[k + 1], β[k + 1] = sum(λ .* x .* p1 .^ 2) / s1, s1 / s
        s = s1
    end
    return α, β
end

"""
    stieltjes_work(w, n, bits; verbose) -> (x, λ, M, gap, usedbits)

The `n`-point Gauss rule for `w`, accurate to `bits`: the discretisation is refined, doubling
the points per piece, until the recurrence coefficients of two successive resolutions agree
to `bits`; then Newton on that recurrence gives the nodes.
"""
function stieltjes_work(w::FunctionWeight, n::Integer, bits::Integer; verbose::Integer = 0,
                        maxpoints::Integer = 4096)
    N = n + 1                                         # the Gauss driver reads a₀ … a_n
    wbits = bits + 32 + 2 * ceil(Int, log2(n + 1))
    tol = ldexp(BigFloat(1), -(bits + 4))
    coefficients(M) = with_bits(wbits) do
        discrete_stieltjes(N, discretize(w, M, wbits)...)
    end
    M = max(2N, 16)
    prev = coefficients(M)
    gap, gaps = BigFloat(NaN), BigFloat[]
    while 2M <= maxpoints
        M *= 2
        cur = coefficients(M)
        gap = with_bits(wbits) do
            scale = max(one(BigFloat), maximum(abs, cur[1]))
            max(maximum(abs, cur[1] .- prev[1]) / scale, maximum(abs.(cur[2] .- prev[2]) ./ cur[2]))
        end
        push!(gaps, gap)
        verbose >= 1 && @info @sprintf("  %d points per piece: recurrence changed by %.2e (want %.2e)",
                                       M, Float64(gap), Float64(tol))
        if gap <= tol
            x, λ, _ = gauss_from_recurrence(n, moment_recurrence(cur...), wbits; name = "StieltjesDiscretization")
            return x, λ, M, gap, wbits
        end
        # A singularity or kink inside g shows as algebraic convergence: the change shrinks by
        # a roughly constant factor per doubling instead of ever faster. After three
        # doublings, extrapolate at the last factor and give up early when the target is out
        # of reach — each doubling costs more than all before it together. A smooth g that
        # starts slowly shrinks by far more than 16 per doubling once M is past 2N.
        if length(gaps) >= 3
            ratio = gaps[end - 1] / gaps[end]
            needed = ratio > 1 ? ceil(Int, Float64(log(gap / tol) / log(ratio))) : typemax(Int)
            if ratio < 16 && (needed > 60 || M * 2.0^needed > maxpoints)
                verbose >= 1 && @info @sprintf("  converging by only %.1f× per doubling: stopping", Float64(ratio))
                break
            end
        end
        prev = cur
    end
    throw(RefinementError("StieltjesDiscretization", """
        the discretisation of $(w) did not converge to $(bits) bits with at most $(maxpoints) points \
        per piece (the recurrence still changed by $(Float64(gap)) at $M points). The discretisation \
        converges quickly only when g is smooth on each piece: move endpoint singularities into the \
        exponents α and β, and split a piece where g has a kink or a jump."""))
end

"""
    StieltjesDiscretization()

Gauss rules for a weight given as a function, a [`FunctionWeight`](@ref): Gautschi's
multiple-component discretisation. Each piece of the weight is discretised with the Gauss rule
of its classical factor, point masses are added as they are, and the recurrence of the
resulting discrete measure comes from the Stieltjes procedure. The number of points per piece
is doubled until two resolutions give the same recurrence to the requested precision, so
the result is accurate as far as the discretisation has converged, and the construction
refuses rather than returns a rule when it does not converge.
"""
struct StieltjesDiscretization <: RuleFamily end

derivation(::Type{StieltjesDiscretization}) = Derived()
family_name(::StieltjesDiscretization) = "StieltjesDiscretization"

const GAUTSCHI_1982 = Citation(key = "Gautschi1982", authors = ["Walter Gautschi"],
                               title = "On generating orthogonal polynomials",
                               journal = "SIAM Journal on Scientific and Statistical Computing", year = 1982,
                               volume = "3", pages = "289--317", doi = "10.1137/0903018")

candidates(::Type{StieltjesDiscretization}, dom::FunctionDomain, ::PolynomialDegree) = [StieltjesDiscretization()]
npoints(::StieltjesDiscretization, dom, degree::Integer) = max(1, cld(degree + 1, 2))
claimed_degree(f::StieltjesDiscretization, dom, degree) = 2 * npoints(f, dom, degree) - 1
degree_range(::StieltjesDiscretization, dom) = 0:typemax(Int)
degree_for_npoints(::StieltjesDiscretization, dom, n::Integer) = 2n - 1
properties(::StieltjesDiscretization, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

function build(f::StieltjesDiscretization, dom::FunctionDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("Gauss nodes for a function weight are irrational; $(T) is not supported"))
    check_support(dom)
    n = npoints(f, dom, degree)
    w = dom.weight
    # k point masses and nothing else: the measure is a k-point rule already, and its
    # orthogonal polynomials stop at degree k, so no n-point Gauss rule with n ≥ k is needed
    # or defined by the recurrence
    isempty(w.pieces) && n >= length(w.atoms) && throw(NoRuleError(
        "$(w) consists of $(length(w.atoms)) point masses only; its Gauss rules have fewer points than that " *
        "(degree ≤ $(2length(w.atoms) - 3)), and the point masses themselves integrate every polynomial " *
        "of degree ≤ $(2length(w.atoms) - 1) exactly"))
    ctx.verbose >= 1 && @info @sprintf("Stieltjes discretisation: %d points for %s, target %d bits", n, w, ctx.bits)
    x, λ, M, gap, used = stieltjes_work(w, n, ctx.bits; verbose = ctx.verbose)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in λ]
    np, na = length(w.pieces), length(w.atoms)
    what = join(filter(!isempty, [np == 0 ? "" : "$np piece$(np == 1 ? "" : "s")",
                                  na == 0 ? "" : "$na point mass$(na == 1 ? "" : "es")"]), " and ")
    cert = Certificate(equations = "recurrence of the discretised measure: $M against $(M ÷ 2) points per piece",
                       residual = BigFloat(gap; precision = 64), residual_bits = used,
                       digits = target_digits(ctx), guard_digits = floor(Int, (used - ctx.bits) * log10(2)))
    prov = Provenance(family = "StieltjesDiscretization", derivation = Derived(),
                      path = ["weight: $w ($what)",
                              "discretisation: Gauss rules of each piece's classical factor, $M points per piece, " *
                              "doubled until the recurrence agreed to $(ctx.bits) bits",
                              "recurrence: discrete Stieltjes procedure at $used bits",
                              "nodes: Newton on the recurrence", "weights: Christoffel function"],
                      seed_source = "none (derived)", citations = [GAUTSCHI_1982])
    return QuadratureRule(xs, ws, dom, PolynomialDegree(2n - 1), prov, cert)
end
