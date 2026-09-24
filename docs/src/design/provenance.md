# Provenance and licensing

Every rule carries a [`Provenance`](@ref CubatureRules.Provenance) record:

```@repl prov
using CubatureRules
p = provenance(rule(Simplex{2}(); degree = 17))
```

| Field | Content |
|---|---|
| `family` | the family that built the rule |
| `derivation` | `Derived()` or `Seeded()` |
| `path` | the steps of the construction: seed, refinement, guard, and any later transformation |
| `seed_source` | where the starting values came from |
| `citations` | the publications the rule is based on |
| `license` | the licence of the numbers in the rule |
| `selection` | why this family was chosen |
| `symmetry` | the symmetry group of the rule, if any |

Operations that change a rule, such as `map_to` or `transform`, add a step to `path`. The
`path` of a rule is therefore a full record of how it was made.

## Citations

`cite(r)` formats the citations of a rule as BibTeX, or as plain text or APA with
`style = :plain` or `style = :apa`. Families attach their citations in `build`. For seeded
families these are the papers that the point counts or the orbit structures come from, even
though the numbers themselves were computed here.

## Licences

The `license` field describes the terms under which the numbers of the rule may be used.
Rules built from the package's own tables, and derived rules, are MIT licensed.

A family whose data comes from elsewhere declares its terms with `family_license(f)`. The
registry stamps this into the provenance of every rule the family builds, in `_build`, so a
family author cannot forget it. A rule refined from such data keeps the licence of the data:
refinement makes the numbers more accurate, but the result is still derived from the
original data.

```julia
julia> provenance(rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 40)).license
"GPL-3.0-or-later, via Lebedev.jl; not covered by CubatureRules.jl's MIT licence"
```

## Selection and licences

A family can also declare `selectable(f) = false`. Its rules are then never returned by
`rule(domain; degree)` unless the caller passes `copyleft = true` or names the family. This
is meant for data with terms, such as the GPL, that a caller should agree to explicitly.
When `rule` passes over a smaller rule for this reason, it prints a warning once for each
case.

## The seed policy

The package does not copy published tables into its own data. All tables in `src/data` were
computed by the package's own searches from orbit structures (see
[Seed strategies](seeds.md)), and `src/data/PROVENANCE.toml` records how each was produced.

Data under other terms can be used through a separate package that provides a family, as
described in [External providers](providers.md). Refining such data to higher precision is
fine, but the result keeps the original licence and is never added to this package's
tables.
