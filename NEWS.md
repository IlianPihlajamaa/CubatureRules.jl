# News

## v0.1.0

The first release. Everything below is new; later entries list changes against the release
before.

### Rules

- `rule(domain; degree)` chooses among every loaded family, fewest points first, and
  `available` and `compare` show the ranking without building anything. `rule(family,
  domain; degree)` names the family.
- Domains: `Interval`, `Simplex` in any dimension, `Orthotope` (boxes), `Sphere`, `Ball`,
  `Disk`, `Wedge`, `Pyramid`, `SurfaceTriangle`, and the unbounded `HalfLine`, `RealLine`,
  `RealSpace`, `GaussianSpace` and `ExponentialSpace`; a `WeightedDomain` pairs a domain
  with a weight (Jacobi, given by moments, logarithmic, point masses, principal values and
  finite parts, oscillatory, and the boundary-element kernels below).
- One-dimensional families: Gauss–Jacobi and its special cases, Gauss–Laguerre,
  Gauss–Hermite, Gauss–Kronrod (through QuadGK.jl), Gauss–Patterson, Fejér, Newton–Cotes,
  tanh-sinh, exp-sinh and sinh-sinh, Filon, and Gauss rules for weights given by moments
  (`ModifiedChebyshev`, `OrdinaryMoments`) or by a discretised measure
  (`StieltjesDiscretization`). Lobatto, Radau and Clenshaw–Curtis come from
  QuadratureRules.jl when it is loaded.
- Fully symmetric rules with positive weights and interior nodes on eleven domains —
  triangle, tetrahedron, 4-simplex, square, cube, 4-cube, disk, ball, wedge, pyramid and
  sphere — stored as `Float64` tables that the package's own searches produced. The point
  counts are the smallest the searches found, not proven minima; on the triangle they are
  below Xiao and Gimbutas's at most degrees from 27 on. No published table is
  redistributed. `LebedevRule` covers every odd degree from 3 to 201 on the sphere.
  `UpstreamLebedev` refines the tables of Lebedev.jl instead; its rules are GPL-3 and are
  only used when asked for.
- Product and combinator families on the other domains: conical products on simplices,
  tensor products on boxes, wedge products, sphere and ball products, Grundmann–Möller.

### Precision

- Any `AbstractFloat` output, `digits = d` for `BigFloat`, and exact `Rational{BigInt}`
  rules where the nodes are rational. A `Float64` request for a tabulated rule is served
  from the table without refinement; any other precision is refined by Newton's method with
  guard digits chosen from the measured condition number.
- The same call gives bitwise identical rules on every platform.

### Checking and records

- `check(r)` verifies exactness against an orthonormal basis of the domain at twice the
  rule's precision, that the rule is not exact one degree higher, and its structure;
  `@test_exact` puts this in a test suite.
- `certificate(r)` records how the refinement went, `provenance(r)` how the rule was made
  and chosen, and `cite(r)` the references. The licence of the underlying data travels with
  the rule.

### Applying rules

- `integrate(f, r)`, and `integrate(f, r, cells)` over a mesh without building a rule per
  cell; `map_to`, `subdivide`, `transform` and `duffy`.
- `EmbeddedRule` and `embedded` for an error estimate at no extra integrand evaluations.

### Boundary elements

- Rules for `1/|y − x₀|`, `1/|y − x₀|³` and its gradient on flat triangles in space, with
  `x₀` on, near or off the panel (`DuffyGauss`, `DuffySinh`, `DuffyFinitePart`); on
  quadratic curved panels (`DuffyCurved`); and in tetrahedra (`DuffyCone`).
- `KernelRuleBuilder` builds the same rules in `Float64` in microseconds, for assembly.

### Known limitations

- Julia 1.11 or later. On 1.11 the `BigFloat` precision is global, so `BigFloat` rules
  should be built from one task at a time; from 1.12 on, it is per task.
