# Makes `LebedevRule{LebedevJLSeeds}` — spelled `UpstreamLebedev()` — constructible from an
# installed Lebedev.jl. Loading that package is all that is needed; the registry finds the
# family with `subtypes` at call time (PLAN §2.5).
#
# Nothing is redistributed here. The numbers come from the caller's own copy of Lebedev.jl,
# and every rule built from them says so: family, licence and citation all point there. The
# selector never takes this variant on its own, and nothing from it is ever written to
# `src/data` — the shipped tables stay in-house.
#
# Refinement is offered. A rule refined from their table is a derived work of it, so the
# GPL-3 licence travels into the refined rule's provenance unchanged. That is the honest
# arrangement: what the seed policy forbids is absorbing someone else's numbers into this
# package's own tables under this package's own licence, not letting a caller refine data
# they already hold.
module LebedevExt

using StaticArrays
# `Lebedev` is also this package's own family name, so the upstream package is imported
# under a name of its own rather than brought in with `using`.
import Lebedev as LebedevPkg
import CubatureRules: build, npoints, claimed_degree, degree_range, cost_estimate, expand,
                      upstream_lebedev_candidates, upstream_lebedev_table, classify_octahedral,
                      refine_octahedral, octahedral_margins, OctahedralMomentSystem,
                      LebedevRule, LebedevJLSeeds, Sphere, PolynomialDegree, isreference,
                      BuildContext, Provenance, Certificate, QuadratureRule, Seeded,
                      NoRuleError, RefinementError, LEBEDEV_LAIKOV_1999,
                      UPSTREAM_LEBEDEV_LICENSE, measure, with_bits, finalize_number,
                      target_digits, isexact, invariant_exponents

# Lebedev.jl indexes its rules by "order", which is the degree of exactness.
_orders() = LebedevPkg.getavailableorders()
_order_for(degree::Integer) = (i = findfirst(>=(degree), _orders()); i === nothing ? nothing : _orders()[i])

# Only the specific method: the generic fallback lives in the main package, and redefining
# it here would be method overwriting, which precompilation forbids.
function upstream_lebedev_candidates(dom::Sphere{3}, c::PolynomialDegree)
    (isreference(dom) && _order_for(c.d) !== nothing) || return LebedevRule[]
    return [LebedevRule{LebedevJLSeeds}()]
end

function upstream_lebedev_table(degree::Integer)
    ord = _order_for(degree)
    ord === nothing && throw(NoRuleError(
        "Lebedev.jl has no rule of degree $degree; it offers $(first(_orders()))–$(last(_orders()))"))
    x, y, z, w = LebedevPkg.lebedev_by_order(ord)
    # Lebedev.jl normalises the weights to sum to 1; here they sum to the domain's measure.
    area = Float64(measure(Sphere{3}()))
    return ord, [SVector{3,Float64}(x[i], y[i], z[i]) for i in eachindex(x)], [area * wi for wi in w]
end

npoints(::LebedevRule{LebedevJLSeeds}, dom::Sphere{3}, degree::Integer) = length(upstream_lebedev_table(degree)[2])
claimed_degree(::LebedevRule{LebedevJLSeeds}, dom::Sphere{3}, degree) = upstream_lebedev_table(degree)[1]
degree_range(::LebedevRule{LebedevJLSeeds}, dom::Sphere{3}) = 0:last(_orders())
degree_range(::LebedevRule{LebedevJLSeeds}, dom) = 1:0
cost_estimate(f::LebedevRule{LebedevJLSeeds}, dom::Sphere{3}, degree, T) =
    T === Float64 ? 1.0 : float(npoints(f, dom, degree))        # a lookup, or a refinement

function build(f::LebedevRule{LebedevJLSeeds}, dom::Sphere{3}, degree::Int, ctx::BuildContext{T};
               seed = nothing) where {T}
    isexact(ctx) && throw(ArgumentError("Lebedev nodes are irrational; $(T) is not supported"))
    ord, xs, ws = upstream_lebedev_table(degree)
    path = ["nodes and weights from the caller's Lebedev.jl, order $ord",
            "weights scaled by the sphere's area: Lebedev.jl normalises them to 1"]
    cert = nothing
    if T !== Float64 || ctx.bits > 53
        # Refine their table with this package's machinery. The orbit structure has to be
        # recovered from the points first, since all that arrives is a cloud.
        got = classify_octahedral(xs, ws)
        got === nothing && throw(RefinementError("Lebedev",
            "the order-$ord table from Lebedev.jl does not decompose into whole O_h orbits, " *
            "so it cannot be refined here; ask for Float64 to get it as it stands"))
        structure, θ64 = got
        θ, res, guard = refine_octahedral(structure, ord, θ64, ctx.bits; cancel = ctx.cancel)
        res.converged || throw(RefinementError("Lebedev",
            "Gauss–Newton did not converge on the order-$ord table from Lebedev.jl " *
            "(residual $(Float64(res.residual)))"))
        wbits = ctx.bits + guard
        bxs, bws = with_bits(() -> expand(structure, θ), wbits)
        xs = [SVector{3,T}(ntuple(j -> finalize_number(ctx, x[j]), 3)) for x in bxs]
        ws = [finalize_number(ctx, w) for w in bws]
        θhat = [finalize_number(ctx, t) for t in θ]
        rbits = 2ctx.bits + guard
        resid = with_bits(rbits) do
            r, _ = OctahedralMomentSystem(structure, ord, BigFloat)(BigFloat.(θhat))
            maximum(abs, r)
        end
        push!(path, "structure recovered from the points: " * string(structure))
        push!(path, "refined here by Gauss–Newton at $(wbits) bits — a derived work of " *
                    "Lebedev.jl's table, and licensed as one")
        cert = Certificate(equations = "O_h-invariant moment system in orbit parameters " *
                                       "(invariants p₄^a p₆^b, 4a + 6b ≤ $ord)",
                           residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                           digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                           cond = res.cond, iterations = res.iterations)
    else
        push!(path, "delivered as it stands; nothing was solved here")
    end
    prov = Provenance(family = "Lebedev", derivation = Seeded(), path = path,
                      seed_source = "Lebedev.jl (installed by the caller, GPL-3)",
                      citations = [LEBEDEV_LAIKOV_1999], symmetry = :Oh,
                      license = UPSTREAM_LEBEDEV_LICENSE)
    return cert === nothing ?
           QuadratureRule(xs, ws, Sphere{3}(), PolynomialDegree(ord), prov) :
           QuadratureRule(xs, ws, Sphere{3}(), PolynomialDegree(ord), prov, cert)
end

end
