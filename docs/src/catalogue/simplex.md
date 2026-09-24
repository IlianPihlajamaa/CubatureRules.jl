# Simplex

```@setup simplex
using CubatureRules, Markdown
function family_table(dom)
    fmt(r) = last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"
    rows = ["| Family | Degrees | Derivation | Symmetry |", "|---|---|---|---|"]
    for row in available(dom)
        der = row.derivation isa CubatureRules.Derived ? "derived" : "seeded"
        push!(rows, "| `$(row.family)` | $(fmt(row.degrees)) | $der | $(row.symmetry) |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

## Domain

```@docs
Simplex
```

On the triangle:

```@example simplex
family_table(Simplex{2}())
```

On the tetrahedron:

```@example simplex
family_table(Simplex{3}())
```

In four and more dimensions:

```@example simplex
family_table(Simplex{4}())
```

## Xiao–Gimbutas

```@docs
XiaoGimbutas
```

These are fully symmetric triangle rules with positive weights and interior nodes, at the
point counts of H. Xiao and Z. Gimbutas, *Comput. Math. Appl.* 59 (2010) 663–676,
doi:10.1016/j.camwa.2009.10.027, or fewer.

The seeds were generated for this package from the orbit structures alone; no published
numbers were used, and the table is MIT licensed. Above degree 20 they were found by growing
from the rule one degree lower and eliminating nodes (see [Seed strategies](../design/seeds.md)),
which found rules with fewer points than the published ones at several degrees, for example
139 instead of 141 at degree 27. Rules with the minimal number of points are not unique; the
table contains the one with the largest smallest barycentric coordinate among those found.

## Fully symmetric tetrahedron rules

```@docs
FullySymmetric
```

Fully symmetric (`S₄`), positive-weight tetrahedron rules with interior nodes. Unlike the
triangle rules, these point counts are not taken from a paper: they are the smallest found by
the package's own search, starting from a number of points below which no fully symmetric
rule can exist. They are not proven to be minimal.

## Grundmann–Möller

```@docs
GrundmannMöller
```

A closed formula for rules of odd degree `2s + 1` on a simplex of any dimension `D`, with

```math
\sum_{i=0}^{s} \binom{s - i + D}{D}
```

points. All nodes and weights are rational, so these are the rules to use for exact
arithmetic. Their weights are negative from degree 3 on, and the number of points grows much
faster than for the other families.

Reference: A. Grundmann and H. M. Möller, *SIAM J. Numer. Anal.* 15 (1978) 282–290,
doi:10.1137/0715019.

## Conical product

```@docs
ConicalProduct
```

The simplex is mapped to a cube by a collapsed (Duffy-type) coordinate transformation, and a
Gauss–Jacobi rule is used in each direction, with Jacobi weights that absorb the Jacobian of
the transformation. With `m = ⌈(d+1)/2⌉` points per direction this gives `m^D` points,
available at every degree and precision.

Reference: A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall
(1971).
