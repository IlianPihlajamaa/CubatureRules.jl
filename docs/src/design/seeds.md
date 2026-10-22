# Seed strategies

A seed is a `Float64` approximation of a rule, accurate enough for Newton's method to
converge from it. Derived families compute their seeds directly; seeded families need a
strategy to find them. The stored seeds are more than that: they are the correctly rounded
parameters of refined rules, from which `Float64` requests are served without refinement
(see [Float64 requests](#Float64-requests)).

## Gauss rules

The Gauss families (Jacobi, Laguerre, Hermite, and rules from moments) take their seeds from
the Golub–Welsch method: the nodes are the eigenvalues of the symmetric tridiagonal Jacobi
matrix of the three-term recurrence, computed in `Float64`. Newton's method on the
recurrence then refines the nodes to working precision, and the weights follow from the
Christoffel function.

## Symmetric rules

The fully symmetric rules are symmetric under a group: on the simplices, the square and the
cube, the disk and the ball, the pyramid, the wedge, and Lebedev's on the sphere. Each is
described by a list of orbits plus a few parameters per orbit (see
[Orbit algebra](symmetry.md)). Its seed is a value for those parameters. All of them are built
the same way, from a row of their seed table. The `seed` keyword of [`rule`](@ref) chooses
where the seed comes from:

| Seed source | Where the starting point comes from | Families |
|---|---|---|
| `TableSeed()` (default) | the table shipped with the package | all |
| `ExplicitSeed(θ)` | orbit parameters you supply | all |
| `MultistartSeed()` | the orbit structure alone, by a multistart search | simplices, sphere |
| `LowerDegreeSeed()` | the family's rule one degree lower, grown and then reduced | simplices |

On the other domains, the searches that produced the tables run from the scripts in
`scripts/`.

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

## Float64 requests

A stored seed is not a rough starting point. Its orbit parameters are those of the rule
refined to high precision, rounded correctly to `Float64`. `scripts/certify_tables.jl` checks
this for every entry: it refines the stored parameters to 160 bits and compares the correctly
rounded result with what is stored, replacing any parameter that is not. (When the first
tables were certified, one of their 78 entries, the degree-11 triangle, was off by 3 ulp.)
The script also records each entry's residual, evaluated at 256 bits, and the test suite
recomputes the recorded residuals and requires them below `10⁻¹⁴`.

So a `Float64` request has nothing left to refine. The stored parameters are expanded into
nodes and weights at 128 bits and rounded once. The weights are then correctly rounded; a
node that is a function of several parameters, `1 − 2a` on the triangle or `√(1 − 2l²)` on
the sphere, is accurate to a few units of rounding, which on a small coordinate can be many
of its own ulps:

```@repl seeds
using CubatureRules
r = rule(Simplex{2}(); degree = 30);
last(provenance(r).path)
@elapsed rule(Simplex{2}(); degree = 30)
```

A `Float32` or `Float16` request rounds the stored rule, and evaluates the residual of the
rounded rule in `Float64`, which resolves it easily. Any request above 53 bits refines the
stored rule as before. Other seed sources (`ExplicitSeed`, `MultistartSeed`,
`LowerDegreeSeed`) are always refined, as is any table entry without a recorded residual.
The test suite recomputes every recorded residual.

## Seeds from other packages

A family can take its seeds from another package, as `UpstreamLebedev` does with
Lebedev.jl. The raw points are decomposed into orbits by
`CubatureRules.classify_octahedral`, and then refined like any other seed. Such rules keep
the licence of their source; see [External providers](providers.md).
