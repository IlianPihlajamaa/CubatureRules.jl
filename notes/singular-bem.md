# Singular and near-singular integration for boundary elements: design note

Working notes for PLAN §0.3 and §6 Tier 5, the flagship application: rules for the singular
and near-singular integrals of boundary-element methods, generated for the element at hand,
at any precision, with a claim that says what they integrate exactly. The spike numbers
below are what the first stage rests on. **Stage 1 is built:**
`WeightedDomain(T, InverseDistance(x₀))` with rules from `DuffyGauss`
(`src/domains/inverse_distance.jl`, `src/families/simplex/duffy_gauss.jl`), as described
under "The core construction", with the API chosen below, on triangles in the plane and in
space (`SurfaceTriangle`, `src/domains/surface_triangle.jl`). A floating-point `x₀` within
rounding of the triangle's plane or boundary is moved onto it: a collocation point computed
in `Float64` is almost never exactly on its element. **Stage 3 (near-singular) is built
too,** as `DuffySinh` (`src/families/simplex/duffy_sinh.jl`), and differs from the plan
below in one respect: the radial direction is not sinh-substituted but gets the Gauss rule
of its own weight on each line, from that weight's moments, which keeps the rule's size
fixed as `x₀` approaches; the sinh substitution is used across the lines only.

## The integrals

A boundary-element code on a triangulated surface needs, for kernels `G(x, y)` that blow up
at `x = y`:

| Case | Integral | Where it comes from |
|---|---|---|
| Collocation, singular | `∫_T φ(y) G(x₀, y) dy`, `x₀ ∈ T` (vertex, edge or interior) | collocation and Nyström methods, the diagonal blocks |
| Collocation, near-singular | the same with `x₀ ∉ T`, at a distance small against the size of `T` | close panels, evaluation near the surface |
| Galerkin, singular | `∫_{T₁} ∫_{T₂} φ(x) ψ(y) G(x, y) dy dx`, `T₁, T₂` equal, or sharing an edge or a vertex | Galerkin assembly |
| Galerkin, near-singular | the same for disjoint but close panels | |
| Volume potentials | `∫_K φ(y) G(x₀, y) dy` on a tetrahedron `K` | volume integral equations |

with `φ, ψ` polynomial (the basis functions on a flat element) and, for Laplace,
`G = 1/|x − y|` (weakly singular on a surface), `∇_y G` (strongly singular: principal value)
and `∂²G` (hypersingular: finite part); Helmholtz multiplies these by `e^{ik|x − y|}`.

## What exists, and what this package should not redo

- **SauterSchwabQuadrature.jl** (in the BEAST.jl ecosystem) implements the
  Sauter–Schwab regularising transformations for the Galerkin cases on triangles and
  quadrilaterals: equal panels, a common edge, a common vertex, positive distance. It takes
  the one-dimensional rule as an argument (a vector of `(x, w)` pairs on `[0, 1]`), so
  arbitrary-precision Gauss rules from this package can be passed straight in. Whether its
  arithmetic is generic enough for `BigFloat` is to be checked, not assumed.
- **Inti.jl** discretises integral operators with
  density interpolation, which regularises the kernel instead of the quadrature, and uses
  Vioreanu–Rokhlin rules on the elements.
- **The literature** for generated rules: the Duffy transform; Johnston and Elliott's sinh
  transformation for near-singular integrals; generalised Gaussian quadrature on Chebyshev
  systems (Ma, Rokhlin and Wandzura 1996) and its nonlinear-optimisation construction
  (Bremer, Gimbutas and Rokhlin 2010), used for the singular integrals of scattering theory
  by Bremer and Gimbutas; Guiggiani's direct method for principal values and finite parts;
  Montanelli, Aussal and Haddar for curved elements.

None of these gives rules generated for a given element at arbitrary precision with an
exactness claim. That is the gap, and it is the one PLAN §0.3 names: for near-singular
integrals the error depends jointly on the kernel and the geometry, so the cases cannot be
tabulated; generation is the only option. This package should provide **rules** — nodes
and weights with a claim and a certificate — and leave assembly, basis functions and
operator algebra to the BEM codes.

## The core construction: a weight in the collapsed direction

Take `x₀` at a vertex `v₀` of a flat triangle `T`, with edges `e₁ = v₁ − v₀`, `e₂ = v₂ − v₀`.
The Duffy map

    y = v₀ + s (e₁ + t (e₂ − e₁)),   (s, t) ∈ [0, 1]²,   dy = s |e₁ × e₂| ds dt

gives `|y − v₀| = s √q(t)` with `q(t) = |e₁ + t (e₂ − e₁)|²`, a quadratic that does not vanish on
`[0, 1]`. So

    p(y) / |y − v₀| dy = p(s, t) |e₁ × e₂| / √q(t) ds dt,

and `p(s, t)` is a polynomial of degree at most `d` in `s` and in `t` when `p` has degree `d`.
The `s` cancels: there is no singularity left, only a weight `1/√q(t)` in the collapsed
direction. Gauss–Legendre in `s` times the Gauss rule of that weight on `[0, 1]` therefore
integrates `p(y)/|y − v₀|` **exactly** for every `p` of degree `≤ d`, with `⌈(d+1)/2⌉²` points.

That is an ordinary `PolynomialDegree(d)` claim on a weighted domain — the triangle with the
weight `1/|y − v₀|` — the case PLAN §6 lists as "product-integration rules for singular
kernels". The weight depends on the triangle's shape through `q`, so the rule has to be
generated per element: one Gauss rule for a weight known in closed form, which the package
already builds at any precision (`FunctionWeight`, or the moments of `t^k/√q(t)`, which
satisfy a three-term recurrence). The other cases reduce to this one:

- `x₀` on an edge or inside `T`: split `T` into two or three triangles with a vertex at
  `x₀`.
- Tetrahedra: the Duffy map from a vertex gives `dy = s² t (…)` and `|y − v₀| = s √q(t, u)`,
  so `1/|y − v₀|` leaves a smooth integrand with a two-dimensional weight `1/√q(t, u)`. That
  weight is not a product, so an exact rule needs a two-dimensional moment problem; until
  then a product rule is spectrally accurate but claims nothing.

### Spike

`notes/singular-bem-spike.jl`. The triangle `(0, 0), (1, 0), (0.3, 0.8)`, singular at the
origin, 40 digits: Gauss–Legendre in
`s` times `rule(WeightedDomain(Interval(0, 1), FunctionWeight(t -> 1/√q(t), 0, 1)); degree = d)`
in `t`, with no new code. The reference is independent of the construction: polar
coordinates about `v₀`, the radial integral `∫₀^{R(θ)} r^{|α|} dr` in closed form and the
angular one by a 200-point Gauss–Legendre rule at 120 digits. Largest relative error over
`y^α/|y|`, `|α| ≤ d`:

| `d` | points | this rule | Duffy × Gauss–Legendre in `t` |
|---|---|---|---|
| 5  | 9   | 9.9e-41 | 5.2e-03 |
| 10 | 36  | 3.7e-40 | 2.1e-05 |
| 20 | 121 | 1.2e-40 | 8.2e-09 |

Exact to the working precision at every degree; for odd `d` it fails at `d + 1` (6.9e-3 at
`d = 5`), and for even `d` the Gauss rules are exact one degree further, as Gauss rules are.
Without the weight, the same points converge only geometrically, because `1/√q(t)` is then
part of the integrand. An ordinary triangle rule is far worse, since the kernel is not a
polynomial: degree 10 (25 points) 1.0e-2, degree 20 (79 points) 4.0e-3, degree 40 (291
points) 8.9e-4.

Helmholtz, `e^{ikr}/r` with `k = 5`, on which this rule is not exact (see below): 4.2e-3 at
9 points, 1.6e-8 at 36, 7.2e-13 at 121, 1.4e-18 at 256 — geometric, from the even-power
terms.

## Beyond the weakly singular vertex case

1. **Helmholtz.** `e^{ikr}/r = Σ (ik)ⁿ rⁿ⁻¹ / n!`: the odd powers are `p/r`, which the rule
   above integrates exactly, but the even ones are polynomials, and the `1/√q` rule in `t` is
   not exact on those. A rule exact on both — `t^b` and `t^b/√q(t)` in the collapsed
   direction — is a generalised Gaussian rule, with a `SpanOf` claim, computed by Newton on
   the moment equations as in Ma–Rokhlin–Wandzura: the package's refinement core applied to
   a non-polynomial system. If the `2n` functions form a Chebyshev system, `n` nodes suffice.
   That has to be checked: a combination `P + Q/√q` vanishes where `P² q = Q²`, a polynomial
   of degree `2n`, one more zero than the Haar condition allows. If it fails, node
   elimination from a larger rule (Bremer–Gimbutas–Rokhlin) still works, with a few more
   nodes.
2. **Near-singular collocation.** With `x₀` at height `h` above the plane of `T` and its
   projection `x₀'` taken as the Duffy apex (with signed sub-triangles when `x₀'` is outside
   `T`), `|y − x₀|² = s² q(t) + h²`. For each `t` the weight `s / √(q(t) s² + h²)` is known in
   closed form, so the `s`-rule can be its Gauss rule, exact in `s`; the result is analytic in
   `t` (since `q > 0`) for every `h ≥ 0`, so Gauss–Legendre in `t` converges geometrically.
   Exact in one direction and empirical in the other: the claim has to say so, or be
   `NoClaim` with a convergence sweep.
3. **Strongly singular and hypersingular kernels.** On a flat panel the double-layer kernel
   vanishes for `x₀` in the plane. The gradient of the single layer gives `1/s` after the
   Duffy map: a finite part in `s` (the package has Hadamard finite-part rules in 1D,
   `FinitePart`), but the finite part with respect to `s` and with respect to the distance
   `r = s √q(t)` differ by a `log √q(t)` term, and the angular parts cancel only when that is
   accounted for (Guiggiani). This needs the most care and comes last.
4. **Galerkin pairs.** Sauter–Schwab reduces each pair configuration to integrals over `[0, 1]⁴`
   with a regular integrand. First check whether SauterSchwabQuadrature.jl runs at `BigFloat`
   with this package's Gauss rules; if it does, document that and test it. Only if it does not
   would this package build its own pair rules — a new domain, a pair of panels.
5. **Curved elements** follow Montanelli–Aussal–Haddar (locate the preimage of `x₀` by Newton,
   then the same constructions in the parameter domain). Later.

## Stages

1. Weakly singular collocation on flat triangles: `x₀` at a vertex, on an edge, inside.
   Rules with a `PolynomialDegree` claim on the kernel-weighted triangle, verified against an
   independent reference (polar coordinates about `x₀`, as in the spike). API to decide.
2. The Helmholtz-ready generalised Gaussian rule in the collapsed direction (`SpanOf`).
3. Near-singular collocation.
4. Galerkin pairs: interop with SauterSchwabQuadrature.jl, or own pair rules.
5. Strongly singular and hypersingular kernels.
6. Tetrahedra; curved elements.

## Open questions

- **API.** A weighted domain such as `WeightedDomain(T, InverseDistance(x₀))`, so that
  `rule(WeightedDomain(T, InverseDistance(x₀)); degree = 10)` returns an ordinary rule whose
  weights carry the kernel, fits the existing design: `integrate(φ, r)` is then
  `∫_T φ(y)/|y − x₀| dy`. The alternative is a dedicated family. The weighted domain is the
  better fit with `PolynomialDegree` claims and with `verify`.
- **Which kernels first.** Laplace `1/r` is the natural first; `r^{-α}` for other `α` is the
  same construction with `q(t)^{-α/2}` and `s^{1-α}` (Gauss–Jacobi in `s`).
- **Scope of the Galerkin case**, as above: delegate if the existing package can be fed
  arbitrary-precision rules.
