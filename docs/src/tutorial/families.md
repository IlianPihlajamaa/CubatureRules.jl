# I want to build a family

A family is a type and a few methods. You don't need to register it: when `rule` or
`available` is called, the package looks up every loaded subtype of `RuleFamily`, including
the ones in your own package.

Here is a complete family with a single one-point rule on any simplex:

```julia
using CubatureRules, StaticArrays
import CubatureRules: candidates, build, npoints, claimed_degree, properties, degree_range,
                      derivation, isreference, Derived, Provenance, NoSymmetry

struct CentroidRule <: RuleFamily end

candidates(::Type{CentroidRule}, dom::Simplex, c::PolynomialDegree) =
    isreference(dom) && c.d <= 1 ? [CentroidRule()] : CentroidRule[]
npoints(::CentroidRule, dom, degree) = 1
claimed_degree(::CentroidRule, dom, degree) = 1      # degree 0 is answered by the degree-1 rule
properties(::CentroidRule, dom, degree) =
    (positive = true, interior = true, symmetry = NoSymmetry(), nested = false)
degree_range(::CentroidRule, dom) = 0:1
derivation(::Type{CentroidRule}) = Derived()

function build(::CentroidRule, dom::Simplex{D}, degree::Int, ctx) where {D}
    T = CubatureRules.outtype(ctx)
    x = SVector{D,T}(ntuple(_ -> one(T) / (D + 1), D))
    w = one(T) / factorial(D)
    prov = Provenance(family = "CentroidRule", derivation = Derived())
    return QuadratureRule([x], [w], Simplex{D}(), PolynomialDegree(1), prov)
end
```

Once your package is loaded, `available(Simplex{2}(); degree = 1)` lists `CentroidRule` and
`rule(Simplex{2}(); degree = 1)` uses it. The package's own test suite includes a separately
compiled package like this one, to make sure this keeps working.

## The methods

| Method | Purpose |
|---|---|
| `candidates(F, dom, claim)` | Returns the instances of `F` that can serve the request, or an empty vector. |
| `npoints(f, dom, degree)` | The number of points, computed without building the rule. Used for ranking. |
| `degree_range(f, dom)` | The degrees the family can produce. |
| `properties(f, dom, degree)` | Whether the weights are positive, the nodes interior, and which symmetry the rule has. |
| `derivation(F)` | `Derived()` if the rule comes from a formula, `Seeded()` if it starts from stored data. |
| `build(f, dom, degree, ctx)` | Constructs the rule. |

`npoints` is called for every candidate before any rule is built, so it should be fast. It
should be exact; a family whose construction decides its own size, refining until two
resolutions agree, returns the fewest points it can have and says so with `exact_npoints`.
`candidates` returns instances rather than the type, so a family with parameters can return
an instance configured for the request, as `GaussJacobi` does with the exponents of a Jacobi
weight.

The other methods are optional; each has a default:

| Method | Default | Define it when |
|---|---|---|
| `claimed_degree(f, dom, degree)` | `degree` | the rule built for `degree` is exact to a higher degree (rounded up to an odd one, or a table entry above) |
| `exact_npoints(f, dom, degree)` | `true` | `npoints` is only the fewest points the rule can have |
| `supports_type(f, T)` | floating-point types | the family can or cannot deliver some number type, such as exact rationals |
| `cost_estimate(f, dom, degree, T)` | points × a precision factor | building costs much more than the point count suggests |
| `weight_amplification(f, dom, degree)` | 1 for positive rules | the weights are signed and the ratio of their absolute sum to their sum is known |
| `degree_for_npoints(f, dom, n)` | none | the family takes a request by number of points |
| `describe_family(f)` | the type name | instances differ in a way the name should show, like `GaussJacobi(1, 2)` |
| `selectable(f)`, `family_license(f)` | `true`, `""` | the rules carry terms this package cannot pass on (see [I want to provide external data](providers.md)) |
| `missing_dependency(f)` | `nothing` | the family needs a package that may not be loaded |
| `needs_degree(f)`, `home_domain(f)` | `true`, `nothing` | the family is requested by something other than a degree, such as a level |
| `inner_families(f)` | none, required for combinators | the family is a `CombinatorFamily` built from others; its terms follow from theirs |


## Writing `build`

The `ctx` argument holds the output type, the target precision in bits, a cancellation token
and a verbosity level. A typical `build` looks like this:

```julia
import CubatureRules: BuildContext, isexact, finalize_number, checkcancel

function build(f::MyFamily, dom, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("MyFamily nodes are irrational; $T is not supported"))
    guard = 24 + 2 * ceil(Int, log2(degree + 1))
    x, w = my_work(degree, ctx.bits + guard)          # compute in BigFloat
    xs = [finalize_number(ctx, xi) for xi in x]       # round to the output type once
    ws = [finalize_number(ctx, wi) for wi in w]
    ...
end
```

Some conventions the rest of the package relies on:

- Compute in `BigFloat` with `ctx.bits` plus some guard bits, and round to the output type
  once at the end with `finalize_number`. Rounding along the way can make results differ
  between platforms.
- If possible, base the number of guard bits on something you measure, such as a condition
  number, and record it in the [`Certificate`](@ref CubatureRules.Certificate).
- In loops that can take a long time, call `checkcancel(ctx.cancel)` so the user can
  interrupt the computation.

## Verification

If your family claims a [`PolynomialDegree`](@ref), [`check`](@ref) needs a basis to verify
it against. Intervals, simplices, boxes, spheres, balls and the weighted domains already
have one. For a new kind of domain, add a `verification_basis` method. If your rules are not
exact for any polynomial space, use [`NoClaim`](@ref) instead of a degree.

## Next steps

- [I want to provide external data](providers.md) if your family uses tables from elsewhere.
- [Selection](../design/selection.md) explains how these methods are used for ranking.
