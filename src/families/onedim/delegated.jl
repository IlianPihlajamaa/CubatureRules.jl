# Families delegated to upstream packages (PLAN §0.1): 1D rules that QuadGK.jl and
# QuadratureRules.jl already produce in generic arithmetic are *not* reimplemented here. They
# are constructed upstream at working precision and converted at the boundary, with the
# provenance recording which package produced the numbers.
#
# QuadGK is a hard dependency: it costs 45 ms to load and gives Kronrod rules, whose
# construction (Laurie's algorithm) there is no reason to rebuild.
#
# QuadratureRules.jl is a *weak* dependency. It pulls in Polynomials and costs a second to
# load, which would double `using CubatureRules` on its own (§7 targets sub-second). Its
# families are therefore defined here but only become constructible when the user loads
# QuadratureRules — at which point the registry finds them with no registration step, which
# is exactly what the `subtypes` design is for.

const QUADRATURERULES_JL = Citation(key = "Kraus_QuadratureRules", authors = ["Michael Kraus"],
                                    title = "QuadratureRules.jl: Julia library collecting quadrature rules",
                                    journal = "Zenodo", year = 2020, doi = "10.5281/zenodo.4310382")
const QUADGK_JL = Citation(key = "Johnson_QuadGK", authors = ["Steven G. Johnson"],
                           title = "QuadGK.jl: Gauss–Kronrod integration in Julia",
                           journal = "GitHub", year = 2013)

"""
    Lobatto()

Gauss–Lobatto rules: both endpoints are nodes, the rest chosen for accuracy. An `n`-point
rule has degree `2n - 3`. Nodes and weights come from
[QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl), so this family is
available once that package is loaded.
"""
struct Lobatto <: RuleFamily end

"""
    Radau(side = :right)

Gauss–Radau rules: one endpoint (`:left` or `:right`) is a node. An `n`-point rule has
degree `2n - 2`. Needs QuadratureRules.jl to be loaded.
"""
struct Radau <: RuleFamily
    side::Symbol
    function Radau(side::Symbol = :right)
        side in (:left, :right) || throw(ArgumentError("Radau side is :left or :right, got :$side"))
        return new(side)
    end
end

"""
    ClenshawCurtis()

Clenshaw–Curtis rules on Chebyshev points of the second kind (endpoints included), with
positive weights at every order. An `n`-point rule has degree `n - 1`, or `n` when `n` is
odd. Needs QuadratureRules.jl to be loaded.
"""
struct ClenshawCurtis <: RuleFamily end

"""
    GaussKronrod()

Gauss–Kronrod rules: the `2n + 1`-point extension of the `n`-point Gauss rule, of degree
`3n + 1` (`3n + 2` for odd `n`), which is what makes an error estimate available at no
extra integrand evaluations. Nodes and weights come from
[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl), which implements Laurie's algorithm in
generic arithmetic.
"""
struct GaussKronrod <: RuleFamily end

derivation(::Type{Lobatto}) = Derived()
derivation(::Type{Radau}) = Derived()
derivation(::Type{ClenshawCurtis}) = Derived()
derivation(::Type{GaussKronrod}) = Derived()
describe_family(f::Radau) = "Radau(:$(f.side))"

for F in (:Lobatto, :Radau, :ClenshawCurtis, :GaussKronrod)
    @eval degree_range(::$F, dom::Interval) = 0:typemax(Int)
    @eval degree_range(::$F, dom) = 1:0
end

# point counts: the smallest n whose degree reaches the request
lobatto_points(d::Integer) = max(2, cld(Int(d) + 3, 2))
radau_points(d::Integer) = max(1, cld(Int(d) + 2, 2))
# Chebyshev points of the second kind need at least two points (they include ±1)
cc_points(d::Integer) = max(2, isodd(d) ? Int(d) : Int(d) + 1)
# the (2n+1)-point Kronrod extension has degree 3n+1 for even n and 3n+2 for odd n
kronrod_degree(n::Integer) = isodd(n) ? 3n + 2 : 3n + 1
function kronrod_halves(d::Integer)
    n = 1
    while kronrod_degree(n) < d
        n += 1
    end
    return n
end

npoints(::Lobatto, dom, degree::Integer) = lobatto_points(degree)
npoints(::Radau, dom, degree::Integer) = radau_points(degree)
npoints(::ClenshawCurtis, dom, degree::Integer) = cc_points(degree)
npoints(::GaussKronrod, dom, degree::Integer) = 2 * kronrod_halves(degree) + 1
claimed_degree(::Lobatto, dom, degree) = 2 * lobatto_points(degree) - 3
claimed_degree(::Radau, dom, degree) = 2 * radau_points(degree) - 2
claimed_degree(::ClenshawCurtis, dom, degree) = (n = cc_points(degree); isodd(n) ? n : n - 1)
claimed_degree(::GaussKronrod, dom, degree) = kronrod_degree(kronrod_halves(degree))

properties(::Lobatto, dom, degree) = (positive = true, interior = false, symmetry = :reflection, nested = false)
properties(::Radau, dom, degree) = (positive = true, interior = false, symmetry = :none, nested = false)
properties(::ClenshawCurtis, dom, degree) = (positive = true, interior = false, symmetry = :reflection, nested = true)
properties(::GaussKronrod, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = true)

candidates(::Type{GaussKronrod}, dom::Interval, ::PolynomialDegree) =
    isreference(dom) ? [GaussKronrod()] : GaussKronrod[]

"""
    upstream_rule(family, n) -> (x, w)

Nodes and weights of the `n`-point rule on `[-1, 1]`, from the upstream package, in BigFloat
at the ambient precision. Defined for the QuadratureRules.jl families by the package
extension.
"""
function upstream_rule end

"QuadGK returns the non-negative half of the symmetric Kronrod rule; mirror it."
function upstream_rule(::GaussKronrod, n)
    xhalf, whalf, _ = QuadGK.kronrod(BigFloat, n)
    x = vcat(xhalf, .-reverse(xhalf[1:(end - 1)]))
    w = vcat(whalf, reverse(whalf[1:(end - 1)]))
    p = sortperm(x)
    return x[p], w[p]
end

upstream_of(::Lobatto) = ("QuadratureRules.jl", QUADRATURERULES_JL)
upstream_of(::Radau) = ("QuadratureRules.jl", QUADRATURERULES_JL)
upstream_of(::ClenshawCurtis) = ("QuadratureRules.jl", QUADRATURERULES_JL)
upstream_of(::GaussKronrod) = ("QuadGK.jl", QUADGK_JL)

delegated_points(f::Lobatto, d) = lobatto_points(d)
delegated_points(f::Radau, d) = radau_points(d)
delegated_points(f::ClenshawCurtis, d) = cc_points(d)
delegated_points(f::GaussKronrod, d) = kronrod_halves(d)

const DelegatedFamily = Union{Lobatto,Radau,ClenshawCurtis,GaussKronrod}

# nothing once the upstream package is loaded and defines `upstream_rule`
missing_dependency(f::DelegatedFamily) = applicable(upstream_rule, f, 1) ? nothing : upstream_of(f)[1]

function build(f::DelegatedFamily, dom::Interval, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("$(describe_family(f)) nodes are irrational; $(T) is not supported"))
    n = delegated_points(f, degree)
    guard = 32 + 4 * ceil(Int, log2(n + 1))
    x, w = with_bits(ctx.bits + guard) do
        upstream_rule(f, n)
    end
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    # the delivered rule is compared with the same rule built at twice the precision
    xhi, whi = with_bits(2ctx.bits + guard) do
        upstream_rule(f, n)
    end
    err = with_bits(2ctx.bits + guard) do
        maximum(max(abs(BigFloat(xs[i]) - xhi[i]), abs(BigFloat(ws[i]) - whi[i]) / abs(whi[i])) for i in eachindex(xs))
    end
    pkg, cite = upstream_of(f)
    cert = Certificate(equations = "nodes and weights from $pkg (agreement with the same rule at twice the precision)",
                       residual = BigFloat(err; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = family_name(f), derivation = Derived(),
                      path = ["$(length(xs)) nodes from $pkg at $(ctx.bits + guard) bits",
                              "converted at the boundary and rounded once to $T"],
                      seed_source = pkg, citations = [cite], symmetry = properties(f, dom, degree).symmetry)
    return QuadratureRule(xs, ws, Interval(), PolynomialDegree(claimed_degree(f, dom, degree)), prov, cert)
end
