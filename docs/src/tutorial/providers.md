# I want to provide external data

Many published tables of quadrature rules have no licence, a copyleft licence, or unclear
terms. This package does not include any of them: all tables in `src/data` were computed for
this package and are MIT licensed.

If you want to make such a table available, put it in a package of your own and define a
family for it. Users who load your package can then use the rules through `rule`, with your
licence attached.

```julia
module MyTables

using CubatureRules
import CubatureRules: candidates, build, npoints, degree_range, properties, derivation,
                      selectable, family_license, Seeded, Provenance

struct MyRules <: CubatureRules.RuleFamily end

selectable(::MyRules) = false
family_license(::MyRules) = "CC-BY-4.0, via MyTables.jl; not covered by " *
                            "CubatureRules.jl's MIT licence"

derivation(::Type{MyRules}) = Seeded()
# candidates, npoints, degree_range, properties and build as for any family

end
```

As with any family, no registration is needed.

## Licence and selection

`family_license` sets the licence text. Every rule that `rule` builds with your family gets
this text in its [`Provenance`](@ref CubatureRules.Provenance), without any further work in
your `build` method:

```julia
julia> provenance(rule(MyRules(), Simplex{2}(); degree = 20)).license
"CC-BY-4.0, via MyTables.jl; not covered by CubatureRules.jl's MIT licence"
```

`selectable(f) = false` means your family is never chosen automatically. It still appears in
[`available`](@ref), and users get it when they name it (`rule(MyRules(), dom; degree)`) or
pass `copyleft = true`. When `rule` picks a larger rule because a smaller one is not
selectable, it prints a warning once for that case.

Use `selectable(f) = false` when your licence has conditions that users may not want to
accept without noticing, such as the GPL. For a permissive licence you can leave it out, and
your rules then compete with the others on number of points.

## Refining the data

Users can ask for your rules at higher precision than your table. The package then refines
your values, and the resulting rule keeps your licence, since it is derived from your data.

The package never stores such refined rules in its own tables. [`UpstreamLebedev`](@ref LebedevRule) is
an example of this setup: it reads the tables of an installed
[Lebedev.jl](https://github.com/stefabat/Lebedev.jl), refines them to the requested
precision, and marks the result as GPL-3.

## Suggestions

- Ship large tables as Julia artifacts rather than files in the repository.
- Record where each table came from and under which terms, as `src/data/PROVENANCE.toml`
  does in this package.
- If the licence of a table is unclear, ask the authors.

## Next steps

- [Provenance and licensing](../design/provenance.md)
- [External providers](../design/providers.md) covers families that delegate computation to
  another package, such as Gauss–Kronrod through QuadGK.jl.
