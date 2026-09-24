# External providers

Some rules are best computed or supplied by another package. There are three ways this
happens.

## Families that delegate computation

One-dimensional rules that other packages already compute in generic arithmetic are not
reimplemented. They are computed by the other package at working precision, converted, and
rounded once like every other rule. The provenance records which package produced the
numbers.

| Family | Computed by | Available |
|---|---|---|
| `GaussKronrod` | [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) | always |
| `Lobatto` | [QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl) | after `using QuadratureRules` |
| `Radau` | QuadratureRules.jl | after `using QuadratureRules` |
| `ClenshawCurtis` | QuadratureRules.jl | after `using QuadratureRules` |

QuadGK loads quickly and is a regular dependency. QuadratureRules.jl pulls in Polynomials.jl
and would roughly double the load time of this package, so it is a weak dependency: the
families are defined here, and a package extension provides their `build` methods once
QuadratureRules is loaded. Before that, `available` still lists the families and `rule`
explains which package to load.

QuadratureRules.jl also exports `nodes` and `weights`. With both packages loaded through
`using`, those names are ambiguous in your session; use `CubatureRules.nodes(r)`, or import
only what you need from one of the two.

## Families that use another package's data

`UpstreamLebedev()` is `LebedevRule{LebedevJLSeeds}()`. When
[Lebedev.jl](https://github.com/stefabat/Lebedev.jl) is loaded, an extension reads its
tables (up to degree 125), decomposes each into octahedral orbits with
`classify_octahedral`, and refines it to the requested precision. Lebedev.jl is GPL-3, so
these rules carry that licence and are not selected automatically; see
[Provenance and licensing](provenance.md).

The extension is activated by loading Lebedev.jl in any way, for example with
`using Lebedev`.

## Families in other packages

Any package can define a subtype of `RuleFamily` and the methods described in
[I want to build a family](../tutorial/families.md). The registry finds it at call time, so
no registration is needed. With `family_license` and `selectable`, such a package can also
supply tables under terms this package cannot include; see
[I want to provide external data](../tutorial/providers.md).

The package's test suite contains a separately precompiled package with its own family, to
make sure such families are found.
