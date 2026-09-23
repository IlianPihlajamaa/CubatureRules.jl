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
is the silent answer to a request that said only "degree 19". `UpstreamLebedev` is the
worked example.

## Publishing a data package

Tables of rules keep turning up in papers whose licence is absent, copyleft, or merely
unstated. This package ships none of them. The supported arrangement is that the data lives
in *your* package, and CubatureRules reaches it through the ordinary family interface:

```julia
struct MyTables <: CubatureRules.RuleFamily end

CubatureRules.selectable(::MyTables) = false          # never chosen on its own
CubatureRules.family_license(::MyTables) = "CC-BY-4.0, via MyTables.jl; not covered by " *
                                           "CubatureRules.jl's MIT licence"
```

plus the usual `candidates`, `npoints`, `degree_range` and `build`. There is no registration
step: the registry finds the family with `subtypes` as soon as your package is loaded.

Declaring `family_license` is all that is needed for the terms to reach the caller — every
rule built through [`rule`](@ref) has them stamped into its `Provenance` automatically. That
is deliberate: leaving it to each author to remember is how a data package ends up producing
rules that *look* unencumbered.

Three things worth doing beyond the interface. Ship large tables as Julia artifacts rather
than repository files — a few tens of megabytes of nodes is a download, not a git history.
Record where the numbers came from and under what permission, as `src/data/PROVENANCE.toml`
does here. And if the licence is unclear, ask the author: permission is cheap, permanent,
and worth more than any workaround.
