# Delegated families

1D rules that other packages already produce in generic arithmetic are not reimplemented
here. They are built upstream at working precision, converted at the boundary, and rounded
once like everything else, with the provenance recording which package produced the numbers.

| Family | Upstream | Availability |
|---|---|---|
| `GaussKronrod` | [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) | always (hard dependency) |
| `Lobatto` | [QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl) | after `using QuadratureRules` |
| `Radau` | QuadratureRules.jl | after `using QuadratureRules` |
| `ClenshawCurtis` | QuadratureRules.jl | after `using QuadratureRules` |

```julia
using CubatureRules, QuadratureRules     # the extension loads with the package

available(Interval(); degree = 7)        # now lists Lobatto, Radau and Clenshaw–Curtis
rule(Lobatto(), Interval(); degree = 11, digits = 50)
```

## Why one is weak and the other is not

QuadGK costs about 45 ms to load, so it is an ordinary dependency, and its Kronrod
construction (Laurie's algorithm) is exactly the kind of thing not worth rebuilding.

QuadratureRules.jl pulls in Polynomials and costs about a second, which would roughly
double `using CubatureRules` against a sub-second target. It is therefore a *weak*
dependency: the families are declared here, and loading QuadratureRules makes them
constructible. Nothing needs registering — the selector finds families with `subtypes` at
call time — so `available` starts listing them the moment the package is loaded. Asking for
one before that gives a message saying which package to load.

## A name clash worth knowing about

QuadratureRules.jl also exports `nodes` and `weights`. Loading both packages unqualified
makes those two names ambiguous in your session:

```julia
using CubatureRules, QuadratureRules
CubatureRules.nodes(r)     # qualify, or import just what you need
```

Everything else in the two packages coexists happily.
