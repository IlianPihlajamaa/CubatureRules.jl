# Adding a rule family

A family is one type plus a handful of methods. There is no registration step: the
selector finds every loaded subtype of `RuleFamily` when it is called.

```julia
using CubatureRules, StaticArrays
import CubatureRules: candidates, build, npoints, properties, degree_range, derivation

struct CentroidRule <: RuleFamily end

candidates(::Type{CentroidRule}, dom::Simplex, c::PolynomialDegree) =
    isreference(dom) && c.d <= 1 ? [CentroidRule()] : CentroidRule[]
npoints(::CentroidRule, dom, degree) = 1
properties(::CentroidRule, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)
degree_range(::CentroidRule, dom) = 0:1
derivation(::Type{CentroidRule}) = Derived()

function build(::CentroidRule, dom::Simplex{D}, degree::Int, ctx) where {D}
    T = CubatureRules.outtype(ctx)
    x = SVector{D,T}(ntuple(_ -> one(T) / (D + 1), D))
    prov = Provenance(family = "CentroidRule", derivation = Derived())
    return QuadratureRule([x], [one(T) / factorial(D)], Simplex{D}(), PolynomialDegree(1), prov)
end
```

After `using` your package, `available(Simplex{2}(); degree = 1)` lists `CentroidRule`
first and `rule(Simplex{2}(); degree = 1)` builds it. The test suite checks exactly this
with a separately precompiled package in `test/downstream`.

A family that ships data also needs an entry in `src/data/PROVENANCE.toml`. A family that
makes an exactness claim needs a `verify` method able to check it. The built-in ones cover
`PolynomialDegree` on intervals and simplices.

A family whose rules carry terms you cannot pass on to your caller — because they come from
a copyleft-licensed source, say — should define

```julia
CubatureRules.selectable(::MyFamily) = false
```

so that `rule(domain; degree)` never returns one on its own. `available` still lists it, and
a caller can have it by naming the family or passing `copyleft = true`; what it will not be
is the silent answer to a request that said only "degree 19". Give `family_license` too, and
record the terms in the rule's `Provenance`, so that they travel with the rule rather than
living in a README. `UpstreamLebedev` is the worked example.
