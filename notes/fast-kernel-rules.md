# Fast kernel rules for boundary-element assembly: design note

`KernelRuleBuilder` (src/fast/) builds the rules of `DuffyGauss`, `DuffySinh` and
`DuffyFinitePart` in Float64 in microseconds, for codes that need one per pair of collocation
point and panel. This note records why it is built the way it is: what was measured, what was
tried and dropped, and how the rules are validated when nothing is verified per rule.

## The problem

A collocation assembly with 10⁴–10⁵ panels needs about 20–50 singular or near-singular rules
per collocation point: 10⁵–10⁶ rules. The certified constructions take 0.03 s (`DuffyGauss`,
degree 15), 0.1–0.35 s (`DuffySinh`), 1–9 s (`DuffyCurved`): hours to days for one assembly,
however good each rule. The cost is the certified path itself — BigFloat moments, Wheeler or
Stieltjes in BigFloat, and an angular refinement that builds the rule two to seven times
until two resolutions agree — which is right for a certified rule at 40 digits and wrong in an
assembly loop for a Float64 rule.

## What is the same and what changes

The constructions are unchanged: the cut at `x₀` (or at the point nearest it), the Duffy map,
the Gauss rule of each ray's or each piece's weight, the finite-part weights folded in along
the rays. Point counts and degree claims are those of the certified families, except that the
near-singular rules are smaller (below). What changes is the arithmetic and the sizing.

**One one-dimensional primitive.** Every weight is, in units of a ray or an edge,
`s^j / ((s − s₀)² + κ²)^{p/2}` on `[0, 1]`: `j = 1` radially in DuffySinh (`p = 1` for `1/r`,
`p = 3` for `1/r³` and the gradient kernel), `j = 0` in the angular direction of DuffyGauss
(`p = 1`) and DuffyFinitePart (`p = 3`). With `s = s₀ + κ sinh z` the weight becomes `dz` or
`dz/(κ cosh z)²`; Gauss–Legendre in `z` is a discretisation, the Stieltjes procedure on it the
recurrence, QL the nodes, a Newton step on the recurrence their refinement, and Christoffel
numbers (a sum of positive terms) the weights. src/fast/quadratic_weight.jl.

**Measured, 2026-10-12** (radial weight, `(s₀, κ)` over `[−10, 0] × [1e-12, 10]`, moment
errors against the certified rule rounded to Float64, the best a Float64 rule can do):

| | worst moment error | µs (m = 5, 10) |
|---|---|---|
| (A) closed-form moments + Chebyshev's algorithm in Double64 or Float64x2 | 1e-14 at m ≤ 5, fails at m ≥ 8 for large `H/q` | 6–15 |
| (A) the same in Float64x3, Float64x4, Float128, BigFloat-128 | x3: 6e-14 to m = 10; Float128, BigFloat fail at m = 10 | 10–74 |
| (B) sinh grid + Stieltjes, all Float64 | 2e-13 – 6e-13 (10–1000× the floor) | 2–8 |
| (B) all Float64x2 | the floor, at every m ≤ 10 | 19–98 |
| (B) Float64, with the grid's `z`-coordinates in double-double | 4e-15 (floor 1e-15) | 2.6–6 |

The arithmetic type barely changes the speed of (A): its cost is allocation and the
eigensolver. (A) fails where the forward moment recurrence loses about 7 bits a step (`x₀`
far from a short ray). (B) in Float64 loses accuracy not in the Stieltjes procedure — every
variant with a Float64 grid was equally bad, whatever precision the later steps used — but in
the grid: `s` grows like `e^z` on an interval up to `|z| ≈ 35`, so a rounding error `eps |z|`
in `z` is a relative error `eps |z|` in `s`, which a moment of degree `k` amplifies `k`-fold;
and the ends must satisfy `κ (sinh z₁ − sinh z₀) = 1`, or every high moment is off by
`(1 − s₀) eps` — 11 ulps for `s₀ = −10`, which was the last factor of ten. With the ends, the
Gauss–Legendre nodes and `z = z_m + h x` in double-double, and `sinh`, `cosh` taking a
first-order correction from the low part, everything else can be Float64. The only extended
operations left are two `asinh` (or `log`) per rule, through MultiFloats.jl: light (SIMD and
MPFR_jll), where DoubleFloats.jl would bring back GenericLinearAlgebra and SpecialFunctions.

**The geometry in Float64x2.** The first validation run found the near-singular `1/r³` wrong
by 1e-7 where the certified rule rounded to Float64 is right to 6e-11: with `x₀` 1e-12
diameters above a triangle whose coordinates are O(1), the height `n·(X − V₁)` computed in
Float64 is wrong by `eps |X|`, a relative error of 1e-4, and `∫ 1/r³ ≈ Ω/h` inherits it. The
certified path computes the geometry exactly in rationals. The fast path now computes the
projection, the nearest point and the offset `w = c − x₀` in Float64x2 and rounds the apex
`c` to Float64: the kernel is evaluated through `w`, consistent with the apex, so a rounded
apex only bends the pieces by an angle of order `eps` and the error is `eps`, not `eps |X|/h`.

## Sizes in advance

**The grid of the primitive.** Calibrated over `s₀ ∈ [−30, 4]`, `κ ∈ [1e-14, 30]` and
`m ≤ 10` for the four weights, against the converged floor (the worst of three large grids):
`M − 8 − 2m` was at most `10.5 L` (`L` half the length of the interval in `z`) for three of
them and `17.4 L` for `q^{-3/2}` in the angular direction; `quadratic_weight_gridsize` takes
`16 + 2m + c min(L, 6) + (c/3) max(L − 6, 0)` with `c = 12` and `20`. The need levels off at
large `L`, since the moments are dominated by the last few units of `z`.

**DuffySinh's angular lines.** After `t = t₁ + κ sinh v`, `t₁ ± iκ` the complex points where
the edge point meets `x₀` (`|a + w + t d|² = 0`; the certified build's `t* ± i√(η² + h²/|d|²)`
is the same point when the nearest point is the foot of the perpendicular, and the edge
version needed fewer lines, 24–28 against 32, beside an obtuse triangle), the angular
integrand is analytic in the strip `|Im v| < π/2`, so Gauss–Legendre converges like `ρ^{-2n}`
with `ρ = b + √(1 + b²)`, `b = π/(2L)`. Against the certified rules (three triangle shapes,
five placements, degree 9) the need was bounded by `4 + K/1.5 + 0.85 ln(1/ε)/(2 ln ρ)`,
`ε = 4e-15`; the prediction alone overestimates long intervals (73 against 56) and
underestimates short ones (17 against 20), where the polynomial degree dominates. The
certified build, refining until two resolutions agree to `2^-(bits + 4)`, delivers two to four
times more lines than Float64 needs.

**The finite part's logarithmic terms** are integrals of `v ln|v| q^{-3/2}` on the same grid,
sized by the larger of the primitive's grid and the angular formula.

## Validation instead of verification

A certified rule carries evidence about itself; a fast rule carries none, so the algorithm is
validated in advance over every geometry it accepts. The rules depend on the geometry only up
to similarity — for a piece, the angle at the apex, the ratio of its sides, and `x₀`'s offset
in units of its size — and scripts/validate_fast_kernel_rules.jl samples that: random
triangles in normalised form from thin to obtuse, rotated, scaled and shifted in space; `x₀`
at vertices, on edges, inside, inside within 1e-7–1e-1 of an edge, above or below at heights
1e-12 to 3 diameters, beside an edge at in-plane distances 1e-6 to 0.3 diameters, and in the
plane outside. For each kernel and degree 0–19, 72 geometries; the fast and the certified rule
(at 20 digits) integrate the orthonormal Dubiner polynomials, and the record
(src/data/fast_kernel_validation.toml) keeps the worst relative error, where it occurred,
and the floor of the certified rule rounded to Float64. `validation(b)` reads it and the
certificate of `rule(b, T, k)` quotes it; `check` still verifies any single rule.

Outside the validated range the builder does not extrapolate: a piece with `η < 1e-8`
(DuffyGauss), `η < 1e-3` (DuffyFinitePart: the `q^{-3/2}` grid's floor degrades as `x₀`
approaches the interior of the far edge, 8e-14 at `η = 1e-4` and 8e-10 at `1e-8`), or a
far-edge distance below `1e-6` of the edge (DuffySinh) goes to the certified construction,
and `b.fallbacks` counts it.

## Results

The record is src/data/fast_kernel_validation.toml. Two findings of the sweep changed the
code before it was recorded: the pieces' edge vectors are now taken from the Float64x2 apex
(with them from the rounded apex, the finite part of `1/r³` at degree 0 was off by 5e-14
where the rounded certified rule is right to 6e-17: the edge vectors moved `x₀` by
`eps |X|` rather than `eps |a|`), and the angular line count has a floor on its degree term
(the gradient kernel at degree 0, 1.3e-12 before, 3.7e-15 with four more lines). For `1/r`
the worst error over 72 geometries is 7e-16 to 2e-14 from degree 0 to 7, within a factor of
three of the rounded certified rule; the finite-part kernels are within that factor at every
degree but 0, where the rounded certified rule is exact to rounding (one node per ray, so
nothing to round but the weights) and the fast rule is backward stable rather than exact: its
weights come from edge vectors rounded once, and the gradient kernel's principal value cancels
between the pieces (1e-13 at worst).

## Not covered

- Curved elements: `DuffyCurved` has no closed-form structure for the rays' weights, and its
  geometry has about fourteen parameters up to similarity, which no sweep covers densely;
  validation within bounds on the curvature would be the way.
- Targets other than Float64: the primitive in Float64x2 reaches the Float64 floor at five to
  ten times the cost, which would be the starting point for Double64 rules.
- Far-field selection: when a plain rule of some degree suffices for `x₀` at a given distance
  is a question for error bounds, not for this builder.
