# CubatureRules.jl — Development Plan

A Julia package providing a large, unified collection of quadrature and cubature rules
across many domains, with an emphasis on **generating** rules on demand at arbitrary
order and arbitrary precision rather than tabulating fixed-precision coefficients.

The name states the scope. One dimension is largely covered by existing packages and is
treated here as a dependency, not a deliverable. What this package owns is dimensions
≥ 2, the symmetry machinery that makes minimal rules possible there, and the
arbitrary-precision generation that makes both tractable.

---

## 0. Ground rules

### 0.1 Name, and the relationship to what already exists

**`QuadratureRules.jl` is not dormant and is not a competitor.** JuliaGNI's package
(Michael Kraus, MIT, DOI `10.5281/zenodo.4310382`) released v0.2.1 on 12 August 2026. It
is 1D only, `BigFloat` by default, and covers Gauss–Legendre, Lobatto, Radau,
Clenshaw–Curtis, Gauss–Chebyshev, tanh-sinh and the tabulated Riemann / midpoint /
trapezoidal rules. Its release notes guarantee that nodes and weights are identical bit
for bit across versions.

That is the interval case, already done, actively maintained, under a compatible licence.
So the policy is not to compete for the name but to divide the work:

- **Do not reimplement** the 1D families it already covers. Depend on it, or convert from
  it at the boundary.
- **Contribute upstream** the 1D families it lacks and that fit its scope — generalised
  Gauss–Jacobi, Laguerre, Hermite, exact-rational Newton–Cotes, Fejér 2, Patterson.
- **`CubatureRules.jl` owns** dimensions ≥ 2, the symmetry and orbit machinery, the
  registry and selector, the arbitrary-measure layer, and the verification subsystem.
- **Open the conversation early.** One note to Kraus before the first release costs
  nothing and avoids two people building the same thing.

The same applies to **QuadGK.jl**, which already ships arbitrary-precision Gauss *and*
Gauss–Kronrod rules for arbitrary weight functions, from a supplied Jacobi matrix or
constructed numerically against the weight. Use `QuadGK.gauss` and `QuadGK.kronrod`
rather than rebuilding Laurie's algorithm. **PolyChaos.jl** is MIT — not GPL, as was
previously assumed here — maintained under SciML, so its Stieltjes and Lanczos
implementations are reusable and its authors are reachable.

Other Julia packages in adjacent territory, none of which does arbitrary precision:
`Lebedev.jl`, `GrundmannMoeller.jl`, `HAdaptiveIntegration.jl`, `Cubature.jl`,
`FastGaussQuadrature.jl`, `GaussQuadrature.jl`. Name them in the documentation's
comparison table; the comparison is favourable, and a registry reviewer will ask.

### 0.2 Licensing and table provenance — the main legal risk

quadpy, the obvious point of comparison, is no longer open source. Current releases on
PyPI carry `LicenseRef-Proprietary` and require a purchased key; the original git history
was deleted from GitHub, and only redistributions of the last GPLv3 release survive
(e.g. `zfergus/legacy-quadpy`). This is precisely the gap this package would fill — but
it also means the most convenient source of tabulated data is GPL'd and cannot be
copied into an MIT-licensed package.

**Policy, enforced from the first commit:**

- Every stored table has an entry in `src/data/PROVENANCE.toml` recording: original
  author, paper DOI, the exact source the numbers were taken from, its license, and the
  date retrieved.
- Prefer authors' own released data (Witherden–Vincent, Xiao–Gimbutas, Womersley,
  Jaśkowiec–Sukumar) and public-domain collections (Burkardt).
- Never transcribe from the GPL quadpy fork, even "just to check".
- CI check: any file in `src/data/` without a provenance entry fails the build.

Retrofitting this later is miserable. It is a small amount of tooling now.

#### Prefer structures to tables

The policy above governs redistributing other people's *numbers*. But the part of a
published rule actually needed as a seed is the **orbit structure** — a small integer
partition such as *"degree 17 on a triangle: one centroid orbit, two vertex-type orbits,
three general orbits"*. That is a combinatorial fact stated in the papers' prose, not an
expressive work, and it is not what a table's licence covers.

Given the structure, multistart Gauss–Newton onto the moment variety recovers the rule
without ever touching the table. Where this works it converts the licensing problem from
a permanent constraint into a compute cost, and it means shipped seed data can be
MIT-licensed from the start rather than only after Tier 4 lands. Test it as part of §0.4;
the machinery is the same either way.

### 0.3 The primary user

**A researcher computing near-singular or high-precision integrals over meshed two- and
three-dimensional geometry** — boundary-element and finite-element work at accuracy that
double precision cannot reach.

Naming one user is not decoration. This plan otherwise has four audiences pulling in
different directions: production FEM kernels want an allocation-free inner loop, people
writing papers want exact rationals and LaTeX, method comparison wants the benchmark
harness, orthogonal-polynomial specialists want the OPQ layer. Each is served, but where
two of them conflict, the BEM/FEM user decides.

The choice is justified by §6 Tier 5: for near-singular integrals the error depends
jointly on the kernel and the element geometry, so the space of cases cannot be
enumerated in a table at all. That is the one place where generation is not an
improvement on tabulation but the only option — where this package's thesis is
structurally rather than merely quantitatively true.

### 0.4 Validate the thesis before building on it

Everything below assumes that Newton on an orbit moment system converges reliably from a
published double-precision seed. That assumption is load-bearing and currently untested.

**The first piece of work is a throwaway script**, not package code: refine exactly one
published rule — Xiao–Gimbutas degree 10 on a triangle — from its `Float64` table to 100
digits, measuring the Jacobian condition number along the way.

If it goes smoothly, the thesis holds and the rest is engineering. If it fights back, the
roadmap needs rewriting before it is committed to. Run the same refinement a second time
from **random** starting points with only the orbit structure held fixed — that tests the
licensing escape route in §0.2 at no extra cost. No other single piece of work in this
plan carries as much information per unit of effort.

---

## 1. Design thesis: seed → refine → certify

Every rule in the package, regardless of family or domain, is produced by the same
three-stage pipeline. Making this the architectural spine — rather than a per-family
grab bag of independent algorithms — is what makes "arbitrary precision, arbitrary
order" tractable, and it is the thing that distinguishes this package from a Julia port
of quadpy.

**Seed.** A cheap `Float64` approximation to the rule. Sources: an asymptotic formula,
a dispatch to `FastGaussQuadrature.jl`, a stored table from the literature, an orbit
structure plus multistart search (§0.2), or a lower-order rule of the same family.

**Refine.** Newton (or Gauss–Newton) on the *defining equations* of the rule, carried
out in generic arithmetic at the requested precision. Quadratic convergence takes 15
digits to 1000 digits in roughly seven iterations. The defining equations differ by
family — roots of an orthogonal polynomial, a symmetric moment system, a nonlinear
system in orbit parameters — but the Newton driver, line search, guard-digit policy
(§2.7), cancellation token and convergence logic are shared.

**Certify.** Evaluate the residual of the defining equations in higher precision (or in
`Arb` ball arithmetic, when available), and attach the measured result to the returned
rule as a `Certificate`.

Note carefully that this last stage is **not** verification (§8). Certifying says *the
Newton iteration converged on the system it was given*. Verifying says *the resulting
rule integrates what it claims to integrate*. A rule seeded from a mis-transcribed table
can have a tiny defining-equation residual and still be the wrong rule — which is exactly
the failure mode §11 worries about. They are different checks, they cost very different
amounts, and they get different types.

This abstraction covers more than it first appears to:

| Rule | Seed | Defining equations |
|---|---|---|
| Gauss–Legendre, order 500, 300 digits | `FastGaussQuadrature` in Float64 | $P_n(x_i) = 0$ via 3-term recurrence |
| Xiao–Gimbutas triangle, degree 20, 200 digits | published Float64 table | symmetric moment system in orbit parameters |
| Xiao–Gimbutas triangle, degree 20, MIT-licensed | published *orbit structure* + multistart | the same system, from random starts |
| Lebedev sphere, order 131, 100 digits | published table | octahedral-orbit moment system |
| Gauss rule for user weight $w$ | discretised Stieltjes in Float64 | orthogonality conditions |

The corollary is an obligation to be honest in the metadata. A rule that only exists
because someone found it by numerical search in 1985 should say so. Publish a clear
taxonomy of **derived** (constructed from first principles at any order) versus
**seeded** (refined from a published table, available only at tabulated orders) and
expose it in the API as `derivation(rule)`.

---

## 2. Types and API

The type hierarchy and the family interface must be settled before many rules are
written. Rules are cheap to add once the shapes are right and expensive to retrofit once
thirty families have baked in a wrong assumption. The first release (§10) exists partly to
force these decisions against a real, hard case rather than settle them in the abstract.

### 2.1 Three orthogonal hierarchies

Domains, rule families, and the resulting rules are separate concepts and get separate
type hierarchies. Conflating them (as many packages do) makes the combinatorics
unmanageable.

```julia
# --- domains -------------------------------------------------------------
abstract type Domain{D,T} end

struct Interval{T}   <: Domain{1,T}; a::T; b::T; end
struct Orthotope{D,T}<: Domain{D,T}; lo::SVector{D,T}; hi::SVector{D,T}; end
struct Simplex{D,T}  <: Domain{D,T}; vertices::SVector{...}; end
struct Sphere{D,T}   <: Domain{D,T}; centre::SVector{D,T}; radius::T; end
struct Ball{D,T}     <: Domain{D,T}; centre::SVector{D,T}; radius::T; end
struct Polytope{D,T} <: Domain{D,T}; ... end        # convex, facet representation
struct Wedge{T}, Pyramid{T}                          # FEM staples

# weighted, possibly unbounded domains — note the concrete base parameter
struct WeightedDomain{D,T,B<:Domain{D,T},W} <: Domain{D,T}
    base::B
    weight::W
end
const HermiteLine   = WeightedDomain(RealLine(), Gaussian())
const LaguerreRay   = WeightedDomain(HalfLine(), Exponential(α))

# --- rule families -------------------------------------------------------
abstract type RuleFamily end

struct GaussJacobi{T}   <: RuleFamily; α::T; β::T; end
struct GaussKronrod     <: RuleFamily end
struct ClenshawCurtis   <: RuleFamily end
struct NewtonCotes      <: RuleFamily; open::Bool; end
struct TanhSinh         <: RuleFamily end
struct GrundmannMöller  <: RuleFamily end
struct XiaoGimbutas     <: RuleFamily end
struct Lebedev          <: RuleFamily end
struct ConicalProduct{F<:RuleFamily} <: RuleFamily; inner::F; end
struct TensorProduct{Fs} <: RuleFamily; families::Fs; end
struct Smolyak{F}       <: RuleFamily; inner::F; end
```

`WeightedDomain`'s base is parameterised concretely rather than typed as the abstract
`Domain{D,T}`; an abstract field there costs a pointer chase on every access for no
benefit.

### 2.2 What a rule claims

```julia
abstract type ExactnessClaim end

struct PolynomialDegree <: ExactnessClaim; d::Int; end   # exact on polys of degree ≤ d
struct SpanOf{B}        <: ExactnessClaim; basis::B; end # exact on a given finite basis
struct NoClaim          <: ExactnessClaim end            # converges, but is exact on nothing
```

`degree::Int` would hardcode the assumption that every rule is characterised by exactness
on a polynomial space. That assumption fails for families scheduled well before v1.0, and
once `rule.degree` is in docstrings and user scripts it cannot be changed.

The forcing case is **tanh-sinh**, whose claim is a convergence rate — $e^{-cn/\log n}$
for functions analytic on the open interval with integrable endpoint singularities — not
an exactness claim at all. With a bare integer field it has to be spelled `-1` or
`typemax(Int)`, and every consumer has to know the sentinel convention.

**The abstraction covers exactness only, and deliberately stops there.** Exactness is a
binary property of a finite-dimensional space. Convergence is an asymptotic property of a
*family* of rules. These are different axes, and a type mixing them will be wrong for any
rule that has one of each. Asymptotic behaviour — the $\omega^{-p}$ order of oscillatory
rules, the double-exponential convergence rate — belongs to the error subsystem (§6,
acknowledged gaps), not here.

**Two apparent cases collapse into `PolynomialDegree` and get no types of their own:**

- *Weighted measures* — log singularities, Jacobi weights, Hermite, and oscillatory
  weights such as $e^{i\omega x}$. The weight belongs to the domain, not the claim; the
  rule is still exact on polynomials against that measure. Filon-type rules land here.
- *Spherical $t$-designs and Lebedev rules.* Spherical harmonics of degree $\le t$ are
  exactly the restrictions of degree-$\le t$ polynomials to the sphere, so these are
  ordinary `PolynomialDegree(t)` rules.

That leaves roughly 90% of the package on `PolynomialDegree`; the abstraction is cheap
precisely because the common case stays simple. Only `SpanOf` is genuinely irreducible,
for systems that are not polynomials against any single weight: generalised Gaussian
quadrature on Chebyshev systems such as $\{1, x, x^2, \log|x|, x\log|x|, \dots\}$.

Two honest caveats about what the hierarchy does *not* buy:

- `NoClaim` carries no information, so verifying such a rule (§8) still needs
  family-specific knowledge of the expected convergence rate. The hierarchy models
  exactness well and models *not*-exactness as a hole that other machinery fills.
- **Tensor products need a decision.** A tensor rule is exact on a tensor-product
  polynomial space, which is larger than the total-degree space a conservative
  `PolynomialDegree` reports. Report the conservative total degree: it is what the
  node-count-sorted selector must compare on, and it is honest. Record the decision so
  `show` does not print something arbitrary.

**Ergonomics are preserved.** `degree(rule)` remains a function, returning the integer for
`PolynomialDegree` claims and erroring usefully otherwise. The `degree = 17` keyword stays
exactly as it is and constructs a `PolynomialDegree` request. Nobody integrating
polynomials on a triangle ever sees the abstraction.

**One case genuinely does not fit and gets a separate type later.** Gauss–Turán and other
rules consuming derivative values change the *functional form* of the rule — nodes carry
an associated derivative order — not its accuracy claim. Give Turán rules their own type
when the time comes rather than widening `QuadratureRule`.

### 2.3 Claim preservation under composition

Every operation that transports or combines rules must state what happens to the claim.
This is the one place the claim abstraction can silently lie, and a silent lie here is
worse than no metadata at all: the rule still returns plausible numbers.

| Operation | Effect on a `PolynomialDegree(d)` claim |
|---|---|
| `map_to(r, dom)` — **affine only** | preserved: degree $d$ |
| `subdivide(r, dom, n)` | preserved: degree $d$ per cell, hence globally |
| `r₁ ⊗ r₂` | conservative total degree $\min(d_1, d_2)$ (§2.2) |
| `transform(r, φ, Jφ)` — general $φ$ | **destroyed** → `NoClaim` |
| `duffy(r)` | **destroyed** → `NoClaim` |

Under a nonlinear $\varphi$ the mapped rule is exact on
$\{p \circ \varphi^{-1}\cdot|J_{\varphi^{-1}}|\}$, which is not a polynomial space. So:

- `map_to` is **affine-only at the type level**. The degree-preserving case must not be
  reachable through the same entry point as the general one.
- `transform` and `duffy` return `NoClaim` unless the caller explicitly asserts a claim,
  in which case the assertion is recorded in `provenance` as an assertion rather than a
  derivation.

Building this in costs almost nothing; adding it later is a breaking change.

### 2.4 Two representations: construction-time and runtime

```julia
struct QuadratureRule{D,T,Dom<:Domain{D},C<:ExactnessClaim,S,W}
    nodes::S                        # Vector{SVector{D,T}} or SVector{N,SVector{D,T}}
    weights::W
    domain::Dom
    exactness::C
    provenance::Provenance
    certificate::Union{Nothing,Certificate}
end

static(r::QuadratureRule)           # → SVector-backed, isbits, provenance in a type param
```

One type cannot serve both a degree-50 `BigFloat` construction and a FEM inner loop
called billions of times. A `Vector` of `SVector` plus a provenance record plus a
certificate is heap-allocated, non-`isbits`, and two pointer chases from the data; for a
6-point triangle rule a hand-written kernel has the nodes in registers.

So storage is a type parameter and `static(rule)` produces the runtime form. The default
`QuadratureRule` is a **construction-time record**: it carries everything needed to
explain, verify and cite itself. The static form is what goes into the hot loop and onto
a GPU. Say which is which in the documentation; the same object cannot honestly be
advertised as both.

A 1D specialisation storing `Vector{T}` rather than `Vector{SVector{1,T}}` is worth the
duplication; 1D is by far the hottest path.

### 2.5 The family interface and the registry

Two functions carry the user-facing promise of this package:

```julia
available(Simplex{2}(); degree = 17)   # what can satisfy this, ranked, without constructing
rule(Simplex{2}(); degree = 17)        # the best candidate, constructed
```

This is the single biggest improvement available over quadpy, which requires the user to
know they want `quadpy.t2.xiao_gimbutas_17()` out of a dictionary of ~1500 opaque names.
It is what makes the package usable by someone who is not a quadrature specialist, and it
is therefore a first-class feature, not a convenience layer bolted on late.

Both are built on one search over the families currently loaded. **Discovery is
reflective; instantiation is explicit.** Keeping those two separate is what makes the
design work.

#### Discovery: `subtypes`

`InteractiveUtils.subtypes(RuleFamily)`, applied recursively through any abstract
intermediate layers, enumerates the families a session has loaded. This is the right
mechanism, and it beats the usual alternative: no global registry vector, no `__init__`
push, no mutable state whose contents depend on load order, and a downstream package that
defines `struct MyFamily <: RuleFamily` is discovered without registering anything.

Two constraints come with it, both easy to respect and both fatal to miss:

- **Never cache the result at module scope.** `const FAMILIES = subtypes(RuleFamily)` is
  evaluated at precompile time and baked into the image, so it silently omits every family
  any downstream package adds. Calling `subtypes` *inside a function* precompiles
  correctly, because the call is not constant-folded — it reads global state at runtime,
  which is exactly what is wanted. The temptation to hoist it for speed is the trap.
- **It is a reflection scan over every loaded module**, so its cost grows with the user's
  session rather than with this package. That is fine once per `rule(...)` call, which is
  the only way it is ever called: §3.1 guarantees construction is a visible subexpression
  the user hoists out of loops. If it ever appears in a profile, cache it behind a cheap
  validity check — never behind a `const`.

Type instability and allocation are **not** constraints here. `rule(...)` returns an
abstractly-typed rule, the search allocates, and neither matters: construction happens
once and is hoisted by design. Paying an allocation and a few microseconds so that users
need not memorise scheme names is an obviously good trade.

#### Instantiation: `candidates`

`subtypes` returns *types*; the search needs *constructible rules*. These differ for three
reasons:

- `GaussJacobi{T}` carries continuous parameters. A type name does not determine an
  instance.
- `ConicalProduct{F}`, `TensorProduct{Fs}` and `Smolyak{F}` are **combinators over other
  families**. `subtypes` returns the combinator and cannot know what belongs inside it.
- One family may offer several distinct candidates on one domain — different orbit
  structures at the same degree, open versus closed variants.

So the candidate set is not `subtypes(RuleFamily)`. It is the closure of the leaf families
under the combinators, restricted to what applies to the domain. `subtypes` supplies the
alphabet; a second method supplies the grammar:

```julia
# each family declares what it can offer on a domain — empty means "not applicable"
candidates(::Type{XiaoGimbutas}, dom::Simplex{2}, c::PolynomialDegree) # → Vector{XiaoGimbutas}
candidates(::Type{GaussJacobi},  dom::Interval,   c::PolynomialDegree) # → Vector{GaussJacobi}

# combinators recurse into the leaves
candidates(::Type{ConicalProduct}, dom::Simplex{D}, c) =
    [ConicalProduct(inner) for F in leaf_families() for inner in candidates(F, Interval(), c)]
```

Returning an empty vector when the family does not apply makes a separate `supports`
predicate redundant — one method instead of two, and no way for the two to disagree.

Everything else a family declares is per-instance metadata, used to rank and to explain:

```julia
npoints(f, dom, degree)          # without constructing the rule
properties(f, dom, degree)       # (positive, interior, symmetry, nested)
derivation(f)                    # Derived() | Seeded()
degree_range(f, dom)
cost_estimate(f, dom, degree, T)
```

#### Ranking, determinism and depth

- Candidates are sorted by node count, ties broken by a documented total order (derived
  before seeded, then family name). **`subtypes`' traversal order is unspecified**, so the
  sort must be explicit or `rule(...)` will not be reproducible across sessions — which
  §5 item 9 requires. The winning choice is recorded in `provenance`, so the selection is
  both reproducible and explicable.
- The combinators make the candidate space infinite in principle
  (`Smolyak{TensorProduct{...}}`). Bound it: combinators are tried to a fixed nesting
  depth, 1 by default, and deeper compositions are reachable only by explicit
  construction.

#### `available` and friends

```julia
available(Simplex{2}())                                        # families × degree range × npoints
available(Simplex{2}(); degree = 17)                           # candidates for one degree, ranked
available(Simplex{2}(); degree = 17, positive = true, interior = true)
compare(Simplex{2}(), 17)                                      # side-by-side properties
```

`available` earns its place from the first release: a triangle at degree 17 already offers
Xiao–Gimbutas, Grundmann–Möller and a conical product, and the choice between them is a
real one — minimal node count versus exact rationals versus guaranteed availability. It
does not need breadth to be useful.

`rule(domain)` with no degree **errors**, listing the candidate families and their degree
ranges. There is no defensible default degree; the right one depends entirely on the
integrand. And note there is no need for a separate `defaultrule(domain)`:
`rule(domain; degree = d)` is already exactly "the rule the selector picks".

When a request cannot be satisfied, say precisely why and what is nearby:

> *No rule with positive weights and interior nodes exists on `Simplex{3}` at degree 22.
> Nearest: XiaoGimbutas degree 21 (78 points), or ConicalProduct degree 22 (1728 points,
> always available).*

### 2.6 Constructing rules

The primary constructor states *requirements*, not a scheme name:

```julia
rule(Simplex{2}(); degree = 17)
rule(Simplex{2}(); degree = 17, T = BigFloat, positive = true, interior = true)
rule(Interval(0, 1); npoints = 40, family = GaussLegendre())
rule(Sphere{3}(); degree = 41, digits = 60)
```

Explicit family selection remains available and is what power users reach for:

```julia
rule(GrundmannMöller(), Simplex{4}(); degree = 11, T = Rational{BigInt})
```

### 2.7 Precision requests and guard digits

Two ways to ask, and both must work:

```julia
rule(dom; degree = 20, T = BigFloat)      # use ambient BigFloat precision
rule(dom; degree = 20, digits = 50)       # pick T and precision automatically
```

The `digits` form is what most users actually want. Internally, every generator runs
inside `setprecision(BigFloat, prec + guard)` and rounds exactly once at the end. Never
depend on ambient global precision inside a generator — pass it explicitly through the
call stack.

**Guard digits are measured, not assumed.** A fixed $O(\log n)$ guard is right for Newton
on a three-term recurrence and wrong for symmetric moment systems. Minimal symmetric rules
sit at or near the degeneracy locus of the moment map — that is *what makes them minimal*
— so the Newton Jacobian is ill-conditioned exactly in the cases this package exists for,
and its condition number grows far faster than $\log(\text{degree})$.

So the refiner estimates the Jacobian condition number during iteration, sets the guard
from it, and **records both the condition number and the guard used in the
`Certificate`**. That turns a silent accuracy failure into a reported number, and gives
users a principled reason to trust or distrust a high-degree rule.

---

## 3. Using rules: the evaluation UX

Construction is half the package. The other half is applying rules pleasantly, and this
is where a Julia package can be markedly nicer than the Python equivalents.

### 3.1 The core call

```julia
integrate(f, rule)                          # rule carries its own domain
integrate(f, rule, domain)                  # map a reference rule onto a domain
integrate(f, rule, mesh)                    # batched over cells, see §3.2
integrate(f, domain; rtol = 1e-30)          # accuracy-targeted, see §3.4
```

**There is deliberately no `integrate(f, domain; degree = 20)`.** Construction is always
visible in the source. Write `integrate(f, rule(domain; degree = 20))` — eight characters
longer, and the rule construction is now a subexpression the reader can see and hoist out
of a loop. Given that the package does no internal caching (§7), the convenience overload's
only real effect would be to hide a repeated, potentially expensive construction inside a
hot loop.

The surviving domain-taking signature is the tolerance-driven one, and the distinction is
principled rather than arbitrary: `rtol` inherently requires building a *sequence* of
rules, so construction is the algorithm rather than an accidental cost.

The documentation's very first example must be the explicit two-line form, so the first
pattern anyone learns is the one that stays correct when they later put it in a loop.

Requirements on the implementation:

- **Allocation-free and inlinable.** The inner loop compiles to what a hand-written FEM
  kernel would produce — which in practice means it runs on `static(rule)` (§2.4).
- **Generic in the integrand's return type.** Complex-, vector-, matrix-, `Unitful`- and
  `StaticArray`-valued integrands work if the accumulation is written generically. Note
  that `sum(w .* f.(x))` allocates; the right primitive is `mapreduce` with an explicit
  zero, and the fiddly part is that the zero is `zero(typeof(w[1]*f(x[1])))`, not
  `zero(T)`.
- **Batched evaluation.** `integrate(f, rule; batch = true)` calls `f` once with the
  full node matrix, for integrands that vectorise or call out to a solver.
- **Threaded and broadcast variants** behind keyword arguments, not separate functions.

### 3.2 Mapping and composition

Reference rules are useless without cheap transport to physical domains.

```julia
r = rule(RefTriangle(); degree = 12)        # constructed once
map_to(r, triangle)                         # affine image, weights scaled by |J|
integrate(f, r, mesh)                       # sum over all cells, rule reused
```

For a mesh with $10^6$ cells this is the difference between a usable package and an
unusable one. The rule is constructed once by the caller, and the per-cell work is an
affine map. Provide the batched form explicitly — do not make users write the loop.

Related, and nearly free once affine mapping exists:

```julia
subdivide(r, domain, n)     # composite rule: n×n subdivision, one rule per piece
r₁ ⊗ r₂                     # tensor product of rules on product domains
transform(r, φ, Jφ)         # general change of variables with user Jacobian
duffy(r)                    # corner-singularity grading on a simplex
```

Each has a defined effect on the exactness claim; see the table in §2.3. The two nonlinear
operations return `NoClaim`.

Composite rules on subdivided domains are the standard answer for integrands that are
smooth but not analytic, and they are cheap to build.

### 3.3 Rule sequences and refinement

Expose families as lazy sequences of increasing accuracy:

```julia
for r in RuleSequence(ClenshawCurtis(), domain)
    @show degree(r), integrate(f, r)
end
```

For nested families (Clenshaw–Curtis, Fejér, Gauss–Patterson, Gauss–Kronrod, sparse
grids) the sequence reuses nodes, so advancing costs only the new evaluations. This
single abstraction buys three things at once: convergence studies, embedded error
estimates, and the accuracy-targeted interface below.

### 3.4 Accuracy-targeted integration

```julia
integrate(f, domain; rtol = 1e-30)
integrate(f, domain; rtol = 1e-30, family = ClenshawCurtis())
```

This is *order*-adaptive, not space-adaptive: it walks a nested rule sequence until the
difference between successive levels meets the tolerance. It is not a replacement for
`QuadGK.jl` or `HCubature.jl` — those are space-adaptive and better for integrands with
localised features.

**Pitch this as high-precision *cubature*.** QuadGK already does adaptive
arbitrary-precision integration on an interval, and does it well. The gap this fills is in
dimensions ≥ 2, where nothing in the Julia ecosystem currently offers order-adaptive
high-precision integration cleanly.

Return an `IntegrationResult` carrying `value`, `error_estimate`, `neval`, and the rule
used, with a `show` method that displays all of it. Never return a bare number from an
adaptive call.

### 3.5 Display

`show(::QuadratureRule)` should print something like:

```
QuadratureRule{2,BigFloat} on Simplex{2}
  family    : XiaoGimbutas (seeded, Newton-refined)
  exactness : polynomial degree 17 (verified sharp)
  points    : 42, all interior
  weights   : all positive, Σw = 0.500000000000000000000000000000
  precision : 100 digits (guard 12, est. cond(J) = 4.7e8)
  residual  : 3.1e-102  (defining equations, 200-digit arithmetic)
  reference : Xiao & Gimbutas (2010), doi:10.1016/j.camwa.2009.10.027
```

This is cheap to implement and disproportionately valuable: it makes correctness
*inspectable* rather than a matter of trust, and it teaches the user the vocabulary of
the package. Note that the residual line names *which* check it reports (§1 versus §8) —
the two must never be conflated in the display.

---

## 4. The benchmark harness

A first-class tool for answering "which rule should I actually use for *my* integrand?"
Together with `available` (§2.5) this is what lets a non-specialist choose well, and it is
a large part of why a user picks this package over a table of numbers.

### 4.1 Interface

```julia
b = benchmark(f, domain;
              rules             = :all,      # or a vector of families
              degrees           = 1:2:31,
              reference         = :auto,     # or an exact value the user supplies
              T                 = Float64,
              max_points        = 10_000,    # per rule
              max_precision     = 256,       # bits, caps BigFloat escalation
              max_time_per_rule = 5.0,       # seconds
              max_total_time    = 120.0)     # seconds, whole run
```

Returns a `BenchmarkResult` — a table with one row per `(family, degree)`:

| family | degree | npoints | t_construct | t_integrate | value | abs_err | rel_err | positive | interior | status |
|---|---|---|---|---|---|---|---|---|---|---|

The table is the product; plotting is a thin layer on top of it. Implement the
`Tables.jl` interface so users can slice, filter and export it themselves.

`status` records *why each family stopped*: `:completed`, `:time_budget`,
`:point_budget`, `:precision_budget`, `:diverged`, `:unsupported`. A curve truncated by a
budget must be drawn as truncated, not left to look as though the family plateaued.

### 4.2 Enforcing the budgets

Julia cannot cleanly interrupt a running computation, so do not try to implement a
timeout by killing work in flight. Two mechanisms instead:

**Escalate and check between degrees.** Run each family in increasing degree, time each
construction, and stop that family when the previous construction exceeded
`max_time_per_rule` — or when a linear fit in log-time predicts the next one will. The
same loop enforces `max_points` and `max_precision`. No interruption is needed and it
covers essentially every realistic case.

**A cancellation token threaded through the Newton driver.** The refine loop is the long
pole and already iterates, so checking a token once per iteration is nearly free and
gives genuine interruptibility for pathological cases — a near-singular moment system that
grinds at one degree. **Build the token into `refine/` in the first release**, even though
the harness itself comes much later: retrofitting it means touching every family.

`max_total_time` is checked between families and truncates the run, reporting which
families were never attempted.

### 4.3 Measurement traps

**Compilation latency.** The first call to each family includes JIT time and will
dominate cheap constructions entirely. Warm up each family at its lowest degree before
timing anything. Where only one sample is affordable, record that in the result rather
than pretending it is a `BenchmarkTools` minimum.

**Timing hygiene.** Since the package does no internal caching (§7), repeated construction
timings are honest by default — one fewer thing to get wrong. Still interpolate
benchmarked values properly (`$`-interpolation under `@belapsed`) and prevent the
integrand call from being optimised away.

### 4.4 The reference value

`reference = :auto` computes ground truth using the package itself: a very high-degree
rule at 2–3× the working precision, cross-checked against a second, structurally
independent family, with agreement to a target number of digits required before the
value is accepted.

If the two disagree, warn loudly and report both. That usually means the integrand is not
smooth — which is itself the most useful thing the benchmark can tell the user, and it
should be surfaced as a headline finding rather than buried in a log message.

This is a satisfying piece of self-consistency: the arbitrary-precision machinery that
justifies the package's existence is also what makes its benchmark tool trustworthy.

### 4.5 Plots

One entry point, one figure, two panels, in a `Makie`/`Plots` package extension — with
the data always available without them:

```julia
plot(b)                    # two panels: error vs npoints | error vs construction time
plot(b; pareto = true)     # same, Pareto frontier highlighted
```

**Left panel: error vs. node count.** The evaluation-cost axis. Since applying a rule is
exactly $\sum_i w_i f(x_i)$, node count is proportional to evaluation time for *any*
integrand — the reader multiplies by their own `time(f)` to get seconds — and it is
machine-independent, so the panel is comparable across people and across runs.

**Right panel: error vs. construction time.** The axis that varies wildly and
unpredictably between families, and the amortizable one: irrelevant if the rule is reused
across $10^6$ mesh cells, decisive if it is built for a single integral. Keeping it in its
own panel lets the reader decide which regime they are in without the plot pretending to
know.

Both panels log–log, one line per family, shared y-axis. Mark negative-weight and
exterior-node rules with a distinct marker — these are disqualifying for many applications
regardless of accuracy. Draw budget-truncated curves with an open terminal marker so a
truncation never reads as convergence.

Resisting a `reuse`-weighted single-axis version of this is deliberate. It would be more
general and less useful: it forces the user to supply a number they usually do not know,
to get a plot they then have to reason backwards from.

### 4.6 Construction-only benchmarking

```julia
benchmark_construction(domain; degrees = 1:2:41, T = BigFloat, digits = 100)
```

Pure rule-generation cost as a function of degree and precision, with no `f` at all. This
falls out of the same harness with the integrand column removed, and it lands **early**:
it is the package's own performance regression test, it answers "are 100-digit degree-40
tetrahedron rules practical?", and it is useful as soon as one domain has two families.
Run it in CI and track it over time.

### 4.7 Caveat, stated in the docs and in the plot subtitle

Benchmark results for one integrand do not generalise. A rule that wins on $e^{-x^2}$
may lose badly on $\sqrt{x}$. The tool is for *choosing a rule for a specific problem*,
not for ranking rules in general. Say this explicitly and prominently; users will
otherwise screenshot the plot and draw the wrong conclusion.

---

## 5. Cheap wins — high value, low implementation cost

Roughly ordered by (value ÷ effort). Ship these opportunistically as the core allows
rather than batching them into one release; several make the early releases substantially
more pleasant to use.

**1. `cite(rule)` → BibTeX.** Provenance is already stored. One function, and academic
users will love it. `cite(rule; style = :bibtex | :apa | :plain)`. Also
`cite(::BenchmarkResult)` to emit every reference used in a comparison.

**2. Node visualisation.** `plot(rule)` with marker size ∝ |weight| and colour by sign,
domain boundary drawn. Instantly reveals negative weights, exterior nodes and broken
symmetry, and it makes excellent README material.

**3. Code emission for FEM kernels.** `emit(rule, :julia)`, `emit(rule, :c)`,
`emit(rule, :cuda)` writing out literal `const` arrays with full-precision decimals.
Lets people bake a 30-digit-accurate rule into a `Float64` kernel with no runtime
dependency on this package at all. FEM and GPU users will adopt the package for this
feature alone, and it is essentially string formatting.

**4. Exact rational output.** `T = Rational{BigInt}` for Newton–Cotes,
Grundmann–Möller, and the low-order symmetric rules whose coefficients are algebraic.
Combined with a `latex(rule)` method, this is the tool people use when writing papers.
Nearly free given generic arithmetic.

**5. `npoints(family, domain, degree)` without construction.** Needed by the selector
anyway; exposing it lets users cost-plan before committing. Similarly `cost_estimate`.

**6. `check(rule)` on demand.** Re-run verification against a user-specified degree and
precision, returning a `Verification`. Useful for rules constructed manually or loaded
from elsewhere, and it makes verification a public feature rather than test-directory
infrastructure.

**7. Embedded error estimates for free.** For nested families, `integrate(f, rule;
error = true)` returns the difference against the next-coarser embedded rule at zero
additional integrand evaluations. Gauss–Kronrod already works this way; generalise the
mechanism to Clenshaw–Curtis, Patterson and sparse grids.

**8. Integrals.jl backend.** Implement the SciML `SciMLBase.solve` interface in an
extension. Makes the package usable inside the whole SciML stack for the cost of one
adapter file.

**9. `hash`, `==`, and deterministic construction.** Rules must be hashable and comparable
so user code can cache them, and `rule(...)` with identical arguments must return
bitwise-identical output across sessions and machines.

This is cheap to *design for* and needs maintaining. With Newton iteration and convergence
tests, `GenericLinearAlgebra` eigensolvers, adaptive guard digits (§2.7), the unspecified
traversal order of `subtypes` (§2.5) and any threading, it requires a dedicated CI job
comparing rule hashes across Julia versions and platforms. `QuadratureRules.jl` makes
exactly this guarantee in its 0.2.1 release notes — that is both the precedent and the
bar. It is also what makes the no-internal-caching decision (§7) painless, since user-side
memoisation becomes a one-liner.

**10. Diagnostic on failure.** When `rule(...)` cannot satisfy a request, say precisely
why and what is nearby (§2.5). Error messages that teach the domain are a real feature.

**11. Unitful support.** Mostly falls out of generic accumulation, but not entirely: node
coordinates carry units and `map_to` scales by a $|J|$ that has units too. Worth a test
and a docs example — as a test, not an assumption.

**12. A `Certificate`-aware `@assert`-style macro** for user test suites:
`@test_exact rule degree` verifying a rule integrates a polynomial basis exactly.

Deliberately *not* filed as cheap despite appearances: moment fitting for cut cells, node
elimination to discover new rules, and generalised Gaussian quadrature for non-polynomial
Chebyshev systems. All are valuable; none is a small job.

---

## 6. Mathematical inventory, tiered by generatability

### Tier 1 — fully derived: any degree, any precision, any dimension

These carry the package's central claim. Make them exhaustive.

| Family | Method | Notes |
|---|---|---|
| Gauss–Jacobi (⊃ Legendre, Chebyshev 1–4, Gegenbauer) | Newton on 3-term recurrence, Float64 seed | dispatch to `FastGaussQuadrature` for `Float64` |
| Gauss–Laguerre / Hermite (generalised) | same | |
| Radau / Lobatto variants of all of the above | Golub's modified Jacobi matrix | one or two prescribed nodes |
| Gauss–Kronrod | **use `QuadGK.kronrod`** | already generic-arithmetic capable; do not rebuild Laurie |
| Clenshaw–Curtis, Fejér 1 & 2 | closed-form weights / generic DCT | nested |
| Newton–Cotes, closed & open | exact integration of the Lagrange basis | exact `Rational{BigInt}` |
| Gauss–Patterson | nested extension | |
| Tanh-sinh, exp-sinh, sinh-sinh | direct formula | essential for `BigFloat` work; `NoClaim` |
| Tensor products | any 1D families | conservative total-degree claim (§2.2) |
| Smolyak sparse grids | combination technique + node deduplication | any 1D family, nested or not |
| Grundmann–Möller (simplex, any $d$, odd degree) | explicit formula | exactly rational; negative weights at high degree |
| Conical / Duffy product (simplex, any $d$) | Gauss–Jacobi in collapsed coordinates | the always-available fallback at any degree |
| Sphere / ball / disk product rules | Gauss–Legendre in $\cos\theta$ × trapezoid in $\varphi$; Jacobi radial weight | exact for spherical harmonics |
| Stroud $E_n^{r^2}$, $E_n^{r}$, $S_n$, $T_n$ | closed form in arbitrary $d$ | |

### Tier 1b — derived, but with boundary nodes

**Convex polytopes** via Lasserre / Chin–Sukumar divergence-theorem reduction sit here,
not in Tier 1 proper. The method gives exact *integrals of polynomials* by facet
recursion; turning that into a point-and-weight rule places nodes **on the facets**, with
counts growing with facet count and dimension.

Such rules are genuinely available at any degree and precision, but the registry's
`interior = true` filter will reject them, and so will any user whose integrand is
singular or undefined on the boundary. Declare `interior = false` in `properties` and say
so in the documentation, rather than promising in §6 something §2.5 will refuse to hand
over. Positive-weight *interior* rules on a general polytope at arbitrary degree require
moment fitting plus node elimination — Tier 4.

### Tier 2 — Gauss rules for arbitrary measures

A modern arbitrary-precision implementation of the parts of Gautschi's OPQ toolbox that
are **not** already available elsewhere:

- **Modified Chebyshev from modified moments** — the headline
- Christoffel modification algorithms: multiply or divide the weight by a linear factor
- Multiple-component discretisation (`mcdis`)
- Singular cases: $\log(1/x)$ weights, algebraic-logarithmic endpoint singularities,
  Cauchy principal value and Hadamard finite-part rules, all reachable via analytically
  known modified moments

**Note what is already covered.** QuadGK ships arbitrary-precision `gauss` and `kronrod`
for arbitrary weight functions, from a supplied Jacobi matrix or constructed numerically
against the weight. PolyChaos.jl (MIT, SciML) implements Stieltjes and Lanczos. Neither
needs rebuilding — depend on them, or contribute to them.

The remaining claim is sharper than the general one it replaces. The moment-based
Chebyshev algorithm is normally dismissed as too ill-conditioned to use, and PolyChaos
explicitly chose Stieltjes and Lanczos over it for that reason. **At 200 digits the
conditioning stops mattering.** That is a concrete, demonstrable payoff from "generate at
any precision" rather than a nicety, and it makes the best headline example in the
documentation.

### Tier 3 — seeded symmetric rules

Xiao–Gimbutas, Witherden–Vincent, Dunavant, Zhang–Cui–Liu, Jaśkowiec–Sukumar for
simplices; Lebedev and Womersley $t$-designs for spheres; Stroud and cyclic rules for
cubes.

Implement an **orbit algebra**: represent a rule as a list of symmetry orbits with a few
free parameters each. On a triangle: the centroid (weight only), a 3-point vertex-type
orbit $(a, w)$, and a 6-point general orbit $(a, b, w)$. Moment equations are written
against a basis of the invariant polynomial ring, and the right-hand sides are *exact
rationals* — the Dirichlet integral formula gives monomial moments over a simplex in
closed form, and monomial moments over $S^{d-1}$ likewise. Newton on the resulting system,
seeded from a published table or from the orbit structure alone (§0.2), recovers the rule
to arbitrary precision.

This is the crux of the package and the least certain part of it. Four decisions it needs,
which §0.4 exists to inform:

- **How invariant bases are computed.** Molien series, an explicit Reynolds operator, or
  hard-coded bases for the handful of groups that actually matter ($S_3$, $S_4$,
  octahedral, icosahedral). The last is probably right, and should be chosen deliberately
  rather than by default.
- **The solve is not square in general.** For minimal rules the Jacobian is rank-deficient
  or near it — that is what minimality means — so use a **rank-revealing solve
  (Gauss–Newton with pseudoinverse, or Levenberg–Marquardt) from the start**, not as a
  fallback for pathological cases.
- **Guard digits from the measured condition number** (§2.7), recorded in the certificate.
- **Feasibility at scale.** Establish early how many free parameters a degree-15 $S_4$
  tetrahedron rule carries, whether its Jacobian is dense or structured, and what one
  refinement costs at 200 digits.

**This is the feature quadpy structurally cannot offer**, and it should lead the README.
A minimal-point Xiao–Gimbutas rule at 200 digits does not exist anywhere else.

### Tier 4 — seed generation

**Node elimination to *discover* new minimal symmetric rules**, rather than merely
refining published ones. The algorithm is a loop over machinery Tier 3 already requires:
take a rule, drop the lowest-weight orbit, Gauss–Newton the remaining orbit parameters
back onto the moment variety, accept if it converges with positive weights and interior
nodes, repeat. The orbit algebra, invariant bases, exact rational moments and
underdetermined Newton solver all have to exist anyway. The hard parts are robustness
tuning and compute time, not new mathematics.

It can wait because it costs *node count*, not capability — conical product rules give any
degree on any simplex, so no user is ever blocked. But three consequences of never
building it are worse than that suggests:

- Coverage stops where the literature stops. Xiao–Gimbutas reaches degree ~50 on
  triangles and ~15 on tetrahedra; beyond that you fall back to product rules with an
  order of magnitude more points.
- Pyramids, wedges, spherical triangles, anisotropic simplex weights and cut cells have
  thin or absent tables.
- **The §0.2 licensing dependency becomes permanent** unless the orbit-structure route
  works.

That last point is also the strongest argument *for* building it: rules generated
in-house ship as MIT-licensed seed tables. An openly-licensed set of minimal symmetric
rules to high degree would be a contribution beyond Julia, and it is what makes
"generate, don't tabulate" fully true rather than half-true.

### Tier 5 — non-polynomial claims

Every family here is a reason the `ExactnessClaim` abstraction (§2.2) exists. Implement
none before v1.0; make sure the type system can describe all of them from the start.

| Family | Claim | Notes |
|---|---|---|
| Generalised Gaussian on Chebyshev systems (Ma–Rokhlin–Wandzura) | `SpanOf` | the only irreducible case; see below |
| Filon, Levin, numerical steepest descent | `PolynomialDegree` on an oscillatory `WeightedDomain` | the $O(\omega^{-p})$ behaviour is error theory, not a claim |
| Exponentially / trigonometrically fitted rules | `SpanOf` | |
| Product-integration rules for singular kernels | `PolynomialDegree` on a weighted domain, usually | `SpanOf` only if the system mixes families |
| Lattice rules, if ever in scope | `SpanOf` | exact on a trigonometric index set, not a degree |
| Gauss–Turán | `PolynomialDegree` | but needs its own *rule* type (§2.2) |

Parameterised weights need no special handling: keep `QuadratureRule` concrete and let
the family carry the parameter, so `rule(Filon(ω), domain; degree = 8)` returns an
ordinary rule. Two practical consequences for the oscillatory families: the weights are
$\omega$-dependent (they are the moments $\int x^k e^{i\omega x}\,dx$, known in closed
form), so a rule cannot be reused across $\omega$; and they are complex, so `T` must admit
`Complex`. Both fall out of generic arithmetic, but the moment computation needs care —
the closed forms suffer catastrophic cancellation at small $\omega$ and need either a
series branch or enough guard digits.

**Generalised Gaussian quadrature for singular kernels is the flagship candidate and
should be promoted ahead of the rest of this tier.** It is what the primary user (§0.3)
actually needs: for BEM, near-singular integrals depend jointly on the kernel and the
element geometry, so the space of cases cannot be enumerated in a table at all.
Generation is not an improvement over tabulation there — it is the only option. It serves
an audience currently hand-rolling this machinery in every BEM codebase, and it reuses the
same Newton-in-high-precision core as everything else. If v1.0 needs one flagship
application to build a paper around, this is a better candidate than degree-50 triangle
rules.

### Acknowledged gaps

Two areas this plan does not address, recorded so they are omissions by decision rather
than oversight:

- **Error theory.** Nothing here exposes a bound or estimate: given derivative or
  analyticity bounds on $f$, what is the actual error for this rule? Peano kernel
  constants are computable with the package's own machinery and would be a genuinely new
  artifact at high order; Bernstein-ellipse bounds are more useful in practice. This
  subsystem also owns everything the claim hierarchy deliberately excludes (§2.2):
  asymptotic orders in $\omega$ for oscillatory rules, and the convergence rate of the
  double-exponential family. Until it exists, those properties are documented in prose
  and not represented in the type system.
- **Non-smooth integrands** get only composite rules and the Duffy transform, despite
  being the case most users actually have.

### Explicitly out of scope

- Adaptive (space-subdividing) integration — that is `QuadGK.jl`, `HCubature.jl`,
  `Integrals.jl`, `HAdaptiveIntegration.jl`. Integrate with them; do not compete.
- QMC and lattice rules — a genuinely different subject. State the boundary in the docs
  and point at `QuasiMonteCarlo.jl`.
- 1D families already covered by `QuadratureRules.jl` and `QuadGK.jl` (§0.1).

---

## 7. Numerical infrastructure

**No caching.** `rule(...)` constructs and returns; it does not memoise. Rules are plain
immutable values, so a user who wants to reuse one binds it to a variable, and a user
with a more complex pattern reaches for `Memoization.jl` or a `Dict` of their own. An
internal cache would buy a modest amount of convenience in exchange for global mutable
state, an eviction policy, precision-dependent keys, thread-safety questions, and
timings that silently lie in benchmarks. Not worth it.

This decision is what makes removing `integrate(f, domain; degree = 20)` (§3.1) the right
call rather than a pedantic one: with no cache, that overload would silently rebuild a
potentially expensive rule on every loop iteration. With the overload gone, no signature
in the package can hide an expensive construction, so the absence of a cache is never a
surprise. The two decisions only work together.

The registry's use of `subtypes` (§2.5) is **not** an exception to this. It reads global
state but writes none, so there is nothing to evict, invalidate or make thread-unsafe —
provided its result is never cached in a `const`.

**Generic linear algebra.** Needed: symmetric tridiagonal eigensolver, LU with partial
pivoting, QR least squares, and a rank-revealing factorisation for the underdetermined
moment systems (§6 Tier 3) — all over arbitrary `Real`. `GenericLinearAlgebra.jl` covers
most of it. Budget for writing an implicit-shift QL routine yourself if its `BigFloat`
path proves too slow around $n \sim 1000$.

**Number types exercised in CI:** `Float32`, `Float64`, `Double64` (`DoubleFloats.jl`),
`BigFloat` at several precisions, `Rational{BigInt}` (exact families only), `ArbFloat`
(certified enclosures).

**Automatic differentiation is a subsystem, not a CI number type.** Differentiable rules
must use the implicit function theorem rather than differentiating through the Newton
iteration, and the IFT derivative depends on the defining equations — so each family needs
a hand-written rule. Budget for it accordingly; `ForwardDiff.Dual` appearing in a test
matrix is not the same thing as supporting it.

**Certified rules** via ball arithmetic are a differentiator nothing in Julia or Python
currently offers: rigorous enclosures of nodes and weights, and therefore a rigorous
error bound on the integral for integrands with known derivative bounds. Behind an
`ArbNumerics` package extension so it is not a hard dependency.

**Load time.** Hard dependencies limited to `StaticArrays`, `GenericLinearAlgebra`,
`LinearAlgebra`, `InteractiveUtils` (for `subtypes`), and whichever of
`QuadratureRules` / `QuadGK` the 1D dependency (§0.1) settles on. Everything else weak.
Target sub-second `using`.

---

## 8. Verification

A first-class subsystem, not a test directory — and distinct from certification (§1).

**Two checks, two types, two names.**

- `Certificate` (§1) records the residual of the **defining equations** at the precision
  used, plus the estimated Jacobian condition number and guard digits. Cheap, produced on
  every construction, and it answers *did Newton converge on the system it was given?*
- `Verification` (this section) records exactness against a **basis of the claimed
  space**. Expensive, produced on demand and in CI, and it answers *is this rule what it
  claims to be?*

The distinction matters because a rule seeded from a mis-transcribed table can have a tiny
defining-equation residual and still be wrong. It matters for cost too: full verification
of a degree-50 rule in three dimensions is roughly 23,400 basis functions × $n$ points at
double the rule's precision. That is a CI campaign, not something every `rule(...)` call
can afford.

**`verify` dispatches on the exactness claim** (§2.2) rather than assuming polynomials.
The methods below are the `PolynomialDegree` implementation, covering ~90% of the package;
`SpanOf` reuses items 1–3 against the supplied basis. `NoClaim` rules have nothing to
verify by exact integration and are checked instead by observed convergence across a node
sweep — the one place verification is empirical rather than certified, which must be
recorded as such in the result.

1. **Exactness.** Every rule tested against a basis of the polynomial space up to its
   claimed degree, residuals computed at 2× the rule's precision. Use *orthogonal* bases
   — Koornwinder–Dubiner on simplices, spherical harmonics on spheres — not monomials.
   Monomial conditioning produces false failures above degree ~15.
2. **Sharpness.** Verify the rule is *not* exact at degree $d+1$, so claimed degrees
   cannot silently drift upward.
3. **Structural invariants.** Weights sum to the domain measure; nodes lie inside the
   domain; the symmetry group is preserved exactly; positivity where claimed.
4. **Cross-validation.** Independently derived rules of the same degree must agree on
   random smooth integrands; agreement with `QuadGK`/`HCubature`; regression against
   archived tables.
5. **Claim preservation.** Every composition operation in §2.3 is tested for the claim it
   reports — including that `transform` and `duffy` report `NoClaim`.

Every family must declare a claim type and have a corresponding `verify` method; a family
whose claim cannot be checked does not ship. Property-based testing over
`(family, degree, T)` triples (`Supposition.jl`) will find more bugs than hand-written
cases, particularly at family and precision boundaries.

**Turn verification failure into an output.** A systematic pass refining every published
table in the literature and verifying the result produces a list of published rules that
do not survive: transcription errors, claimed degrees that are wrong, rules asserted to be
minimal that are not. That artifact does not currently exist, it costs nothing beyond what
Tier 3 already does, and it is a better paper hook than "we wrote a package" — as well as
the most persuasive available demonstration that generating beats tabulating. §11 files
refinement failure as a risk; it is also a deliverable.

---

## 9. Repository layout

```
src/
  domains/          # geometry, affine maps, measures, exact monomial moments
  symmetry/         # orbit types, invariant bases, group actions
  refine/           # generic Newton, Gauss-Newton, continuation, guard digits,
                    #   cancellation token, seed sources
  families/
    onedim/         # jacobi + whatever is not delegated upstream
    measures/       # OPQ layer: modified chebyshev, christoffel, mcdis
    simplex/  sphere/  cube/  unbounded/  polytope/
  composition/      # tensor, smolyak, conical, duffy, subdivision, polytope reduction
  registry/         # candidates, selection, ranking, available(), compare()
  apply/            # integrate, mapping, batching, static rules, rule sequences
  verify/           # exactness, sharpness, certificates, verifications
  emit/             # code generation, LaTeX, BibTeX
  benchmark/        # harness, reference values, result tables
  data/             # seeds, orbit structures + PROVENANCE.toml
ext/
  MakieExt.jl  PlotsExt.jl  ArbNumericsExt.jl  IntegralsExt.jl
  MeshesExt.jl  SymbolicsExt.jl  UnitfulExt.jl  ForwardDiffExt.jl
```

---

## 10. Roadmap

Ordered by dependency and by value, not by calendar. The principle is that the first
public release should be the thing that cannot be obtained anywhere else, so the hardest
machinery gets built while abandoning the project is still cheap — and so the type
hierarchy and family interface are forced against a real hard case rather than designed in
the abstract.

**v0.1 — the cornerstone, and the interfaces it forces.**
Settle the type hierarchy and the family interface first: domains, `ExactnessClaim`,
`QuadratureRule` with its storage parameter and `static`, `Certificate` and
`Verification`, claim preservation under composition (§2.3), and the
`candidates` / `npoints` / `properties` / `derivation` interface (§2.5).
Then the simplex machinery that proves the thesis: exact rational monomial moments, the
orbit algebra and invariant bases, the Newton/Gauss–Newton refiner with rank-revealing
solve, measured guard digits and the cancellation token.
**Triangle rules at arbitrary precision** — Xiao–Gimbutas refined, Grundmann–Möller exact
rational, conical product as the always-available fallback — plus Gauss–Jacobi in 1D,
because the conical product needs it.
`integrate`, affine `map_to`, `show`, `cite`, `available`, `compare`, the selector,
claim-dispatched verification, `benchmark_construction`. Documentation site live.

*One-sentence pitch: a minimal-point Xiao–Gimbutas degree-20 triangle rule at 200 digits,
which exists nowhere else.*

**v0.2 — tetrahedra and meshes.**
Extend the orbit algebra to three-dimensional simplices, with the invariant-basis and
conditioning questions of §6 Tier 3 answered at $S_4$. Mesh integration at scale, the
`static(rule)` path, batched and threaded application. This is the release the FEM
audience adopts.

**v0.3 — one dimension, completed and delegated.**
Take the dependency on `QuadratureRules.jl` and `QuadGK.jl` decided in §0.1, contribute
upstream what fits there, and fill the remainder here: Laguerre, Hermite, Radau/Lobatto,
Clenshaw–Curtis, Fejér, exact-rational Newton–Cotes, Patterson, tanh-sinh. Tensor
products and boxes. The registry gains real breadth and `available` becomes genuinely
discriminating.

**v0.4 — spheres, balls, unbounded domains.**
Product rules, Lebedev with refinement, $t$-designs, Stroud $E_n^{r^2}$ families.

**v0.5 — arbitrary measures.**
The OPQ layer in arbitrary precision, scoped to what QuadGK and PolyChaos do not already
cover (§6 Tier 2): modified Chebyshev from modified moments, Christoffel modification,
`mcdis`, and the singular-weight catalogue. Oscillatory weights.

**v0.6 — composition and seed generation.**
Sparse grids and Smolyak; wedges, pyramids, general polytopes via divergence-theorem
reduction (Tier 1b); composite and subdivided rules. Begin node elimination (Tier 4),
targeting an MIT-licensed seed table for at least one domain where the literature runs
out — which both extends coverage and starts unwinding the §0.2 licensing dependency.

**v0.7 — ecosystem and tooling.**
The full benchmark harness with integrands and plots; `Integrals.jl`, `Meshes.jl` and
Makie extensions; `ArbNumerics` certified rules; code emission; AD rules.

**v1.0.**
Stable API, ≥95% coverage, complete provenance manifest, JOSS or arXiv paper. Freeze the
query interface; leave rule families open for extension.

Ship the cheap wins from §5 opportunistically throughout rather than batching them into
one release — several of them (`show`, `cite`, diagnostics, `npoints`) make the earlier
releases substantially more pleasant to use.

---

## 11. Risks

**The orbit algebra may be harder than assumed.** It is the crux of the thesis and the
least analysed part of this plan. §0.4 exists to find out before anything is built on it.
If Newton on the moment system does not converge reliably from published seeds, the
roadmap above is wrong from v0.1 onward and needs rewriting, not patching.

**Table provenance** is the largest legal exposure. Build the licence manifest on day one,
and test the orbit-structure route (§0.2) early — it is the difference between a permanent
dependency on other people's terms and a compute cost.

**Newton refinement will not always converge.** Some tabulated rules sit near singular
points of the moment map; some published tables contain typos; some orbit structures give
underdetermined systems. The documented fallback is Gauss–Newton with pseudoinverse or
Levenberg–Marquardt — which for minimal rules is the expected path, not the exception
(§6 Tier 3). Flag tables that fail to refine in the registry; never silently return the
unrefined seed.

**`BigFloat` performance will disappoint** above $n \sim 500$ unless allocation is
controlled and the generic eigensolver path is avoided where a Newton-on-recurrence path
exists. `benchmark_construction` lands in v0.1 precisely so this is discovered early
rather than at v0.5.

**Scope creep.** The Tier-4 and Tier-5 items are seductive and each could absorb the whole
project. Ship Tiers 1–3 completely first.

**Bus factor.** Given quadpy's trajectory, a single-maintainer package is a real risk to
users who would be building on it. Register under an organisation — JuliaMath or
JuliaApproximation are both plausible homes — rather than a personal account, and write a
CONTRIBUTING guide that makes adding a new rule family a well-trodden path. The
`candidates` interface (§2.5) exists partly to make that path short: a new family is one
type, one `candidates` method and its metadata, with no registration step and nothing
central to edit.
