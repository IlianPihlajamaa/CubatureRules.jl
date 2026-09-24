# I want to cite a rule

Each rule records the publications it is based on. [`cite`](@ref) returns them as BibTeX:

```@repl cite
using CubatureRules
r = rule(Simplex{2}(); degree = 17);
print(cite(r))
```

Other formats are available through `style`:

```@repl cite
print(cite(r; style = :plain))
print(cite(r; style = :apa))
```

## How the rule was made

[`provenance`](@ref) returns the full record:

```@repl cite
provenance(r)
```

It lists where the starting values came from, how they were refined and at what precision,
and why this family was chosen over the others. This is useful when you need to describe
exactly which rule you used, for example in a paper.

## Licence

The provenance also contains the licence of the numbers in the rule:

```@repl cite
provenance(r).license
```

Rules built from this package's own tables are MIT licensed. A rule built from data supplied
by another package carries that package's licence instead. For example, a Lebedev rule
refined from the tables in Lebedev.jl is GPL-3:

```julia
julia> provenance(rule(UpstreamLebedev(), Sphere{3}(); degree = 29)).license
"GPL-3.0-or-later, via Lebedev.jl; not covered by CubatureRules.jl's MIT licence"
```

This matters if you publish the nodes and weights themselves. Such rules are never chosen
automatically: you only get them if you ask for them by name or pass `copyleft = true`. See
[Provenance and licensing](../design/provenance.md).

## Next steps

- [I want to provide external data](providers.md)
- [Provenance and licensing](../design/provenance.md)
