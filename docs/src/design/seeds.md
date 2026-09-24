# Seed strategies

A seed is a `Float64` approximation of a rule, accurate enough for Newton's method to
converge from it. Derived families compute their seeds directly; seeded families need a
strategy to find them.

## Gauss rules

The Gauss families (Jacobi, Laguerre, Hermite, and rules from moments) take their seeds from
the Golub–Welsch method: the nodes are the eigenvalues of the symmetric tridiagonal Jacobi
matrix of the three-term recurrence, computed in `Float64`. Newton's method on the
recurrence then refines the nodes to working precision, and the weights follow from the
Christoffel function.

## Symmetric rules

The minimal rules on the triangle, tetrahedron and sphere are symmetric under a group, and
are described by a list of orbits plus a few parameters per orbit (see
[Orbit algebra](symmetry.md)). Their seed is a value for those parameters. The `seed` keyword
of [`rule`](@ref) chooses where it comes from:

| Seed source | Where the starting point comes from |
|---|---|
| `TableSeed()` (default) | the table shipped with the package |
| `ExplicitSeed(θ)` | orbit parameters you supply |
| `MultistartSeed()` | the orbit structure alone, by a multistart search |
| `LowerDegreeSeed()` | the family's rule one degree lower, grown and then reduced |

### Multistart search

Given only the orbit structure, random parameters are drawn and refined with
Levenberg–Marquardt in `Float64`. Starts that converge to a rule with positive weights and
interior nodes are kept. This is enough for low degrees, but becomes slow as the number of
parameters grows.

### Growing from the rule below

`LowerDegreeSeed` starts from the rule one degree lower, which already integrates all but
the highest moments exactly. It adds orbits, refits the parameters, and then removes points
again: it drops an orbit, or merges two parameters so that an orbit becomes a smaller type,
and refits after each move. Several such chains run with different random choices, and the
rule with the fewest points wins.

The shipped tetrahedron seeds above degree 10, the triangle seeds above degree 20 and the
Lebedev seeds above degree 15 were found this way. At several triangle degrees the search
found rules with fewer points than the published tables.

## Where the tables come from

All tables in `src/data` were generated for this package by the searches above, starting
from orbit structures only. No published tables were copied. The point counts are
sometimes taken from the literature as a target, and the literature is cited, but the
numbers are the package's own and are MIT licensed. `src/data/PROVENANCE.toml` records how
each table was produced.

The scripts in `scripts/` regenerate the tables.

## Seeds from other packages

A family can take its seeds from another package, as `UpstreamLebedev` does with
Lebedev.jl. The raw points are decomposed into orbits by
`CubatureRules.classify_octahedral`, and then refined like any other seed. Such rules keep
the licence of their source; see [External providers](providers.md).
