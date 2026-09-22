# Exp-sinh and sinh-sinh: the double-exponential family on unbounded domains
# (PLAN §2.2, §6 Tier 1), the siblings of tanh-sinh in `tanhsinh.jl`.
#
#   exp-sinh   ∫₀^∞ f(x) dx   x = exp(π/2 · sinh t),   w = h · π/2 · cosh(t) · x
#   sinh-sinh  ∫_ℝ f(x) dx    x = sinh(π/2 · sinh t),  w = h · π/2 · cosh(t) · cosh(π/2 · sinh t)
#
# Both report `NoClaim` for the same reason tanh-sinh does — they are exact on no
# polynomial space — and both are nested in the level, so consecutive levels give an error
# estimate at no extra evaluations.
#
# Truncation is the one real design choice. On a finite interval tanh-sinh can stop where
# the weight falls below the resolution of the output precision, because the missing tail
# is then provably negligible. Here the weights *grow* double-exponentially with t, so no
# such bound exists without knowing how fast f decays; the transformation's premise is that
# it decays. The rules therefore truncate at a fixed dynamic range instead: the nodes span
# [2^-2(p-2), 2^2(p-2)] in magnitude, where p is the output precision in bits.
#
# The *square* of the resolvable range, not the range itself, because what the truncation
# drops is a head and a tail, not a weight. An integrand behaving like x^(-1+δ) at the
# origin leaves behind x_min^δ, and one decaying like x^(-1-δ) leaves x_max^(-δ); the worst
# case worth serving is δ = 1/2, and squaring the range puts that at roundoff — ∫₀^∞
# x^(-1/2) e^(-x) dx goes from 5e-8 to 1e-15 in Float64, for 18% more nodes, because t
# enters through asinh. Still, these rules are valid only for integrands that decay, and
# the quantity to watch is the convergence sweep, not the certificate.

"""
    ExpSinh(level = 4)

Exp-sinh (double exponential) quadrature on the half line `[0, ∞)`, with step `h = 2^-level`
in the transformed variable. Suited to integrands that decay at infinity and may have an
integrable singularity at the origin, `∫₀^∞ x^{-1/2} e^{-x} dx` being the model case.

Nodes span `[2^-2(p-2), 2^2(p-2)]` in magnitude at output precision `p`; the head and tail
beyond are assumed negligible, which is exactly the assumption that `f` decays. The range is
chosen so that an integrand behaving like `x^(-1/2)` at either end leaves a remainder at
roundoff. Reports [`NoClaim`](@ref), so it is never
offered for a `degree` request and must be asked for by name,
`rule(ExpSinh(5), HalfLine())`. Verify it with `CubatureRules.verify_convergence` over a
sequence of levels; consecutive levels are nested, so [`EmbeddedRule`](@ref) gives an error
estimate for free.
"""
struct ExpSinh <: RuleFamily
    level::Int
    function ExpSinh(level::Integer = 4)
        1 <= level <= 20 || throw(ArgumentError("exp-sinh level must be between 1 and 20, got $level"))
        return new(Int(level))
    end
end

"""
    SinhSinh(level = 4)

Sinh-sinh (double exponential) quadrature on the whole line, with step `h = 2^-level` in the
transformed variable. Suited to integrands that decay in both directions, `∫_ℝ e^{-x²} dx`
being the model case.

The same truncation caveat as [`ExpSinh`](@ref) applies: nodes span `[2^-2(p-2), 2^2(p-2)]`
in magnitude and the tails beyond are assumed negligible. Reports [`NoClaim`](@ref); ask for it
by name, `rule(SinhSinh(5), RealLine())`.
"""
struct SinhSinh <: RuleFamily
    level::Int
    function SinhSinh(level::Integer = 4)
        1 <= level <= 20 || throw(ArgumentError("sinh-sinh level must be between 1 and 20, got $level"))
        return new(Int(level))
    end
end

const DoubleExponentialUnbounded = Union{ExpSinh,SinhSinh}

derivation(::Type{ExpSinh}) = Derived()
derivation(::Type{SinhSinh}) = Derived()
describe_family(f::ExpSinh) = "ExpSinh(level $(f.level))"
describe_family(f::SinhSinh) = "SinhSinh(level $(f.level))"
family_name(::ExpSinh) = "ExpSinh"
family_name(::SinhSinh) = "SinhSinh"
needs_degree(::DoubleExponentialUnbounded) = false
degree_range(::DoubleExponentialUnbounded, dom) = 0:-1                 # no polynomial degree at all
claimed_degree(::DoubleExponentialUnbounded, dom, degree) = -1
properties(::ExpSinh, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = true)
properties(::SinhSinh, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = true)
npoints(f::DoubleExponentialUnbounded, dom, degree = 0) = length(de_work(f, 53 + 32, 53)[1])

home_domain(::ExpSinh) = HalfLine()
home_domain(::SinhSinh) = RealLine()

# Never answer a `PolynomialDegree` request with a rule that is exact on nothing.
candidates(::Type{ExpSinh}, dom::HalfLine, ::PolynomialDegree) = ExpSinh[]
candidates(::Type{SinhSinh}, dom::RealLine, ::PolynomialDegree) = SinhSinh[]

"""
    de_work(family, bits, outbits) -> (x, w)

Exp-sinh or sinh-sinh nodes and weights in BigFloat of precision `bits`, truncated to the
dynamic range that `outbits` can resolve. Nodes ascend; `t = k h` runs symmetrically about
zero, so halving `h` retains every node of the previous level.
"""
function de_work(f::DoubleExponentialUnbounded, bits::Integer, outbits::Integer = bits)
    return with_bits(bits) do
        h = ldexp(BigFloat(1), -f.level)
        halfπ = BigFloat(π) / 2
        # |x| ≤ 2^2(outbits-2) for both maps, since sinh(u) ≈ e^u/2 for large u
        tmax = asinh(2 * (outbits - 2) * log(BigFloat(2)) / halfπ)
        K = floor(Int, tmax / h)
        K > 1 << 22 && throw(RefinementError(family_name(f), "the node sum did not terminate"))
        xs = Vector{BigFloat}(undef, 2K + 1)
        ws = Vector{BigFloat}(undef, 2K + 1)
        for (i, k) in enumerate(-K:K)
            t = k * h
            u = halfπ * sinh(t)
            c = h * halfπ * cosh(t)
            if f isa ExpSinh
                x = exp(u)
                xs[i], ws[i] = x, c * x
            else
                xs[i], ws[i] = sinh(u), c * cosh(u)
            end
        end
        xs, ws
    end
end

"""
    de_residual(f, x̂, ŵ, bits, outbits) -> BigFloat

How far the emitted numbers are from the closed form they claim to be: the largest relative
deviation of the rounded nodes and weights from `x = exp(π/2 sinh kh)` (or `sinh`) and its
weight formula, recomputed at precision `bits`. `outbits` must be the precision the rule was
emitted at, so that the same nodes are compared.
"""
function de_residual(f::DoubleExponentialUnbounded, xhat, what, bits::Integer, outbits::Integer)
    x, w = de_work(f, bits, outbits)
    length(x) == length(xhat) || throw(ArgumentError("the reference sum has a different length"))
    return with_bits(bits) do
        res = zero(BigFloat)
        for i in eachindex(x)
            res = max(res, abs(BigFloat(xhat[i]) - x[i]) / abs(x[i]),
                      abs(BigFloat(what[i]) - w[i]) / w[i])
        end
        res
    end
end

function build(f::DoubleExponentialUnbounded, dom::Domain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("$(family_name(f)) nodes are irrational; $(T) is not supported"))
    guard = 32 + 2f.level
    x, w = de_work(f, ctx.bits + guard, ctx.bits)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    span = floor(Int, 2 * (ctx.bits - 2) * log10(2))
    cert = Certificate(equations = f isa ExpSinh ?
                                   "x = exp(π/2 sinh kh), w = h π/2 cosh(kh) x (largest relative deviation)" :
                                   "x = sinh(π/2 sinh kh), w = h π/2 cosh(kh) cosh(π/2 sinh kh) " *
                                   "(largest relative deviation)",
                       residual = BigFloat(de_residual(f, xs, ws, ctx.bits + guard, ctx.bits); precision = 64),
                       residual_bits = ctx.bits + guard, digits = target_digits(ctx),
                       guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = family_name(f), derivation = Derived(),
                      path = ["h = 2^-$(f.level), $(length(xs)) nodes",
                              "truncated at 1e-$(span) ≤ |x| ≤ 1e$(span), set by the output precision",
                              "claim: none — the rule converges but is exact on no polynomial space"],
                      seed_source = "none (closed form)", citations = [TAKAHASI_MORI_1974],
                      symmetry = properties(f, dom, degree).symmetry)
    return QuadratureRule(xs, ws, home_domain(f), NoClaim(), prov, cert)
end
