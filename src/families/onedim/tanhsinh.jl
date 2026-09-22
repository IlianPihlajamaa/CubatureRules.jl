# Tanh-sinh (double exponential) quadrature (PLAN §2.2, §6 Tier 1).
#
#     x = tanh(π/2 · sinh(t)),   w = h · (π/2) · cosh(t) / cosh²(π/2 · sinh(t)),  t = k h
#
# This is the family that forces the `ExactnessClaim` hierarchy to exist. Tanh-sinh is not
# exact on any polynomial space: its claim is a *convergence rate*, e^{-cn/log n} for
# functions analytic on the open interval with integrable endpoint singularities. It
# therefore reports `NoClaim`, is not offered for a `PolynomialDegree` request, and is
# verified by a convergence sweep rather than by exact integration (§8).
#
# The nodes crowd towards ±1 double-exponentially, which is what makes the endpoint
# singularities harmless: the weights vanish faster than the integrand blows up. Halving h
# keeps every previous node, so consecutive levels are nested and give an error estimate at
# no extra evaluations.
#
# Upstream note: QuadratureRules.jl also has tanh-sinh, but its interface takes a level whose
# nodes saturate at ±1 in the precision it is given, which does not fit the precision policy
# here (§2.7), so this one is built locally from the closed form.

"""
    TanhSinh(level = 4)

Tanh-sinh (double exponential) quadrature on an interval, with step `h = 2^-level`. Nodes
crowd towards the endpoints double-exponentially, which makes integrable endpoint
singularities harmless — the case Gauss rules handle badly.

Accuracy is limited by how well `1 - x` is represented: the outermost node sits a few units
of roundoff from the endpoint, so an integrand that is evaluated at `x` near ±1 — say
`1/sqrt(1-x^2)` — loses about half the working digits (3e-8 in `Float64`, 3e-26 at 50
digits). Ask for more digits, or substitute so the integrand is written in terms of the
distance to the endpoint.

Reports [`NoClaim`](@ref): it is exact on no polynomial space, so it is never offered for a
`degree` request and must be asked for by name, `rule(TanhSinh(5), Interval())`. Verify it
with `CubatureRules.verify_convergence` over a sequence of levels. Consecutive levels are
nested, so [`EmbeddedRule`](@ref) gives an error estimate for free.
"""
struct TanhSinh <: RuleFamily
    level::Int
    function TanhSinh(level::Integer = 4)
        1 <= level <= 20 || throw(ArgumentError("tanh-sinh level must be between 1 and 20, got $level"))
        return new(Int(level))
    end
end

derivation(::Type{TanhSinh}) = Derived()
describe_family(f::TanhSinh) = "TanhSinh(level $(f.level))"
needs_degree(::TanhSinh) = false
home_domain(::TanhSinh) = Interval()
degree_range(::TanhSinh, dom) = 0:-1                      # no polynomial degree at all
properties(::TanhSinh, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = true)
npoints(f::TanhSinh, dom, degree = 0) = length(tanh_sinh_work(f.level, 53 + 32, 53)[1])
claimed_degree(::TanhSinh, dom, degree) = -1

const TAKAHASI_MORI_1974 = Citation(key = "TakahasiMori1974", authors = ["Hidetosi Takahasi", "Masatake Mori"],
                                    title = "Double exponential formulas for numerical integration",
                                    journal = "Publications of the RIMS, Kyoto University", year = 1974,
                                    volume = "9", pages = "721--741", doi = "10.2977/prims/1195192451")

# A `PolynomialDegree` request is never answered by a rule that is exact on nothing.
candidates(::Type{TanhSinh}, dom::Interval, ::PolynomialDegree) = TanhSinh[]

"""
    tanh_sinh_work(level, bits) -> (x, w)

Tanh-sinh nodes and weights on `[-1, 1]` at `h = 2^-level`, in BigFloat of precision `bits`.
The sum is truncated where both the weight and the distance to the endpoint fall below the
resolution of that precision — beyond which further terms cannot change the result.
"""
function tanh_sinh_work(level::Integer, bits::Integer, outbits::Integer = bits)
    return with_bits(bits) do
        h = ldexp(BigFloat(1), -level)
        # Truncate at the *output* precision: a node any closer to ±1 would round to ±1
        # there, and an integrand with an endpoint singularity would then be evaluated on
        # the singularity itself.
        cutoff = ldexp(BigFloat(1), -(outbits - 2))
        halfπ = BigFloat(π) / 2
        xs = BigFloat[]
        ws = BigFloat[]
        k = 0
        while true
            t = k * h
            u = halfπ * sinh(t)
            c = cosh(u)
            w = h * halfπ * cosh(t) / c^2
            x = tanh(u)
            (w < cutoff && 1 - abs(x) < cutoff) && break
            if k == 0
                push!(xs, x); push!(ws, w)
            else
                push!(xs, x); push!(ws, w)
                pushfirst!(xs, -x); pushfirst!(ws, w)
            end
            k += 1
            k > 1 << 20 && throw(RefinementError("TanhSinh", "the node sum did not terminate"))
        end
        xs, ws
    end
end

function build(f::TanhSinh, dom::Interval, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("tanh-sinh nodes are irrational; $(T) is not supported"))
    guard = 32 + 2f.level
    x, w = tanh_sinh_work(f.level, ctx.bits + guard, ctx.bits)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    # the truncated sum of the weights differs from the measure by the tail that was dropped
    tail = with_bits(2ctx.bits + guard) do
        abs(sum(BigFloat.(ws)) - 2)
    end
    cert = Certificate(equations = "closed-form double-exponential nodes (deviation of Σw from the measure, " *
                                   "i.e. the truncated tail)",
                       residual = BigFloat(tail; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = "TanhSinh", derivation = Derived(),
                      path = ["h = 2^-$(f.level), $(length(xs)) nodes after truncation",
                              "claim: none — tanh-sinh converges but is exact on no polynomial space"],
                      seed_source = "none (closed form)", citations = [TAKAHASI_MORI_1974], symmetry = :reflection)
    return QuadratureRule(xs, ws, Interval(), NoClaim(), prov, cert)
end
