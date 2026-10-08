# Pooled search for fully symmetric rules

A long-running search for new and smaller seeded rules for every fully symmetric family: the
triangle, tetrahedron and 4-simplex; the square, cube and 4-cube; the disk and ball; the pyramid
and wedge; and the sphere. Many single-threaded workers share one pool of exact rules on disk, instead of each
running its own chains of eliminations. This follows the procedure that took the octahedral
sphere rules from degree 191 to 201, generalized to any domain with an adapter.

It never writes to `src/data`. Finds are refined to 100 digits and collected in
`campaign/pooled/<domain>/`. They reach the shipped tables only after
`scripts/certify_tables.jl` has checked them (see [Merging finds](#merging-finds)).

## Running it

```sh
# start workers; they are independent processes and keep running after the launcher exits
julia --project scripts/search_campaign/launch.jl triangle:extend:2 triangle:improve:1 \
      tetrahedron:extend:2 square:extend:1 cube:extend:1

# what has been found, and the state of the pool
julia --project scripts/search_campaign/status.jl

# stop every worker after its current task
touch campaign/pooled/STOP          # or campaign/pooled/<domain>/STOP for one domain
```

A DOMAIN may be a comma-separated list, `disk,ball:extend:1`: that worker takes the domains in
turn, an hour each (`--slice H`), which keeps seven workers enough for eleven domains.

To run a single worker, for testing or by hand:

```sh
julia --project --heap-size-hint=2G scripts/search_campaign/campaign.jl triangle --worker 1 --mode extend
```

**Modes.**
- `extend` works one degree step above the best table, which is the shipped table merged with
  the finds so far. First, though, it settles the frontier: it goes back to the lowest new
  degree whose best rule is still more than three unknowns above square (typically a freshly
  grown rule), until that rule is pruned or the worker has spent two hours there without
  improving it. Otherwise every later degree would be grown from a poor base: in a first run
  that moved up as soon as any rule existed, the triangle reached degree 69 in a few hours,
  but with 673 points at degree 58 against 576 at degree 57.
- `improve` works at degrees between `--from` and `--to` (default: the top four of the best
  table), looking for rules with fewer points.
- `fill` closes a gap in a table, working upwards from `--from` to `--to`: at the lowest
  missing degree (grown from the one below) or the lowest unsettled one, as extend does; once
  the gap is closed, it improves within it. It is what the sphere uses for degrees 19–131 (see
  [The sphere](#the-sphere)).
- `--blind N` makes a worker ignore the shipped entries from degree `N` up. That's how the
  benchmarks below time a rediscovery.

**Memory.** Each worker needs about 1–2 GB; `--heap-size-hint` stops Julia's garbage collector
from letting it grow further. Before the hint was set, seven workers on a 32 GB machine,
alongside other work, ran out of memory.

## The procedure

**Background.** All rules here are fully symmetric: a list of orbits (an orbit structure),
each with a weight and a few coordinates, constrained by the invariant moment equations of
degree `n`. That gives `m` equations for `u` unknowns.

- A rule with `u − m = e > 0` ("excess `e`") lies on an `e`-dimensional family of rules.
- A rule with `e = 0` is *square* and isolated.
- To remove points, you change the structure (drop an orbit, merge two of its values, or send
  one to zero), which removes unknowns, and then refit.

**What the old search does.** The package's own search (`grow_and_eliminate`,
`box_grow_and_eliminate`) runs independent *chains*. Each chain:
1. grows the rule one degree step down by enough orbits to be well above square (about `m/4`
   unknowns), then
2. removes points one move at a time, refitting after every candidate move, until no move
   succeeds.

**Where its time goes.** Profiling triangle degrees 30 and 40 and tetrahedron degree 14:
- 98% of the time went into failed refits, most of them in the last steps of a chain. There,
  almost every move fails, and every candidate is tried before the chain gives up.
- Every chain repeats that endgame from scratch, and chains share nothing.
- About 10% of the failed refits had in fact converged, only onto a rule with a negative weight
  or a node outside the domain.

The pooled search keeps every intermediate rule, so the endgame becomes shared work. A worker
repeatedly picks one of five tasks:

1. **Chain** (25% of the time). One of the package's own chains, from the best rule one degree step
   down. Where a chain costs seconds (boxes, lower degrees) these are hard to beat; see the
   benchmarks.
2. **Grow.** From the best rule one degree step down, add orbits until the structure is some
   number of unknowns above square:
   - the target is chosen at random from 3, 6, 12, `m/8` and `m/4`;
   - in 30% of grows, one orbit of every smaller type is added as well.

   The new orbits are placed at random with small weights, then fitted. A fitted rule goes into
   the pool.
3. **Prune.** From a pool rule more than three unknowns above square, remove points as
   `eliminate` does: the moves that remove most points first, in random order among those,
   refitting after each, as far as any move succeeds. Every intermediate rule goes into the
   pool. The rule where pruning stops is marked *exhausted*: the remaining moves from it are
   left to the finish step.
4. **Finish.** The endgame, shared between workers and retried with variations:
   - **Pick a parent:** a rule at most three unknowns above square, or an exhausted one. Attempts
     are spread evenly over orbit families (the counts of each orbit type). Within a family,
     the rules with the fewest points and the fewest attempts so far come first.
   - **Pick a move:** 70% of the time one that removes the most points.
   - **Perturb:** shift the coordinates by a scale cycling through 0, 1e-4, 5e-4, 2e-3, 6e-3 and
     1.5e-2.
   - **Fit:** alternately with the package's fit and with the positive-weight solve below.
   - A success that still has unknowns to spare goes back into the pool.
5. **Kick** (basin hopping). Prune and finish only ever remove points, so a rule they have
   exhausted is a dead end. A kick leaves it instead:
   - each worker keeps a current rule per degree, at first (and in 10% of kicks again) the best
     one known;
   - it adds one to three orbits (70% general ones) at random places, fits, and prunes back as
     far as any move succeeds, which often removes a different orbit than the one added;
   - the result becomes the current rule if it has no more points, and otherwise with
     probability exp(−ΔN / T), with T half a general orbit's worth of points (3 on the triangle,
     12 on the tetrahedron, 4 on the square, 24 on the cube).

   This is a Metropolis walk over local minima: it can cross a plateau of equal point counts, or
   climb out of a minimum, where every other task only descends.

**Choosing a task.** Apart from the 25% of chains, a worker chooses by what the pool holds at
its target degree:
- *grow* while the pool is nearly empty;
- *prune* while nothing there is near square;
- otherwise a *kick* a quarter of the time, and mostly *finish* the rest, with some grow.

**Candidates.** Any degree-`n` rule with fewer points than the best known, at any excess, is a
candidate:
- it is refined to 100 digits by the package (`refine_seed`);
- it is checked to be positive and interior, with no two nodes closer than 1e-6 (as
  `scripts/verify_tables.jl` requires; closer nodes mean two orbits are merging, which is a
  move for the prune and finish tasks rather than a rule);
- it is written to `candidates/`, and the merged best table is rewritten.

**Positive-weight solve.** Levenberg–Marquardt with every weight written as `w̄ exp(x)`, so that
no weight can change sign along the way:
- steps in the log-weights are capped at 2;
- the solve stops once 40 accepted steps have improved the residual by less than 2%;
- a converged result is handed to the package's own fit for polishing and acceptance, so a rule
  is accepted by exactly the same test either way.

**The pool on disk.** Each rule is stored as `pool/degNNN/eXXX_pYYYYYY_hHASH.toml`:
- The name carries the excess, the point count, and a hash of the structure plus the parameters
  rounded to 1e-7. A rule found twice (which happens often, across several workers) is stored
  once, and a worker can choose tasks from file names alone.
- Files are written under a temporary name and then renamed, so no worker ever reads half a
  file.
- Nothing needs a lock: workers only add files.

## Benchmarks

One pooled worker against one process running the package's chains back to back, each for
15 minutes on one core. The pooled worker ran blind at the target degree, so it had to rediscover
the rule from the degree below, just as the chains did. "Table" is the shipped point count.

| Degree | Table | Chains: best (time to reach it) | Pooled: best (time to reach it) |
|---|---|---|---|
| triangle 35 | 226 | 226 (84 s), 6 chains in total | **225** (7 min; it had reached 226 at about 90 s) |
| tetrahedron 16 | 250 | 263 (4.5 min), 6 chains | **261** (3 min) |
| square 25 | 120 | 120 (43 s), 119 chains | 121 (2 min)¹ |
| cube 13 | 160 | **154** (2 min), 269 chains | **154** (3 min) |

¹ This run predates the chain task. With it, a pooled worker reaches 120 points in 1.5 minutes.

**What the benchmarks show.**
- Pooling pays off where a single chain is expensive: triangles and tetrahedra, and the higher
  degrees generally. At the triangle and tetrahedron frontier (degrees 54–56 and 23–27) a single
  chain costs 10–20 minutes.
- The cheap box chains remain worth running; hence the chain task.
- Both methods beat the cube rule then shipped at degree 13: 154 points against 160. The
  table has had a 154-point rule since.

### Kicks, and two ideas that did not help

Three additions were tried as options against the engine without them, on the same 15-minute,
single-core, blind runs, with different seeds (worker numbers):

- **kick**: basin hopping, as described above;
- **adaptive**: tasks chosen by Thompson sampling of each task's rate of results within one
  general orbit of the best, per second;
- **ucb**: the finish task choosing orbit families by UCB1 on refit success, instead of evenly.

| Problem (table) | Without | kick | kick + ucb | ucb | adaptive |
|---|---|---|---|---|---|
| tetrahedron 16 (250) | 250, 265, 258, 250 | **248**, 259, **248**, 259 | –, –, 259, 259 | 250, 259 | 250, 265 |
| cube 19 (380) | 384, 386 | 374, 376 | 378, **369** | | |
| triangle 40 (291) | 292, 292 | 292, 291 | 291, 291 | | |
| square 41 (305) | 304, 304 | 304, 304 | 305, 304 | | |
| triangle 35 (225) | 226, 226 | 226, 226 | | 226, 226 | 226, 226 |

- **Kicks help.** On cube 19 they end about 10 points lower, on the tetrahedron 2–7 points
  lower on average, and on triangle 40 they reach the table more often. They are never worse.
  They are now a regular task.
- **ucb** is mixed: the best cube result, but the tetrahedron runs were no better than kicks alone.
- **adaptive** changed nothing: in these runs the best counts came from chains and grows in the
  first minutes, before there was anything to adapt to.

Neither of the last two was kept. Three of these runs beat the shipped tables (tetrahedron 16:
248 points against 250; cube 19: 369 against 380; square 41: 304 against 305).

## Layout

```
campaign/pooled/
  STOP                        optional: stops every worker
  workers.txt                 workers started by launch.jl, with process ids
  <domain>/
    pool/degNNN/*.toml        exact rules above square, Float64
    candidates/*.toml         finds, refined to 100 digits, one table row each
    <table>.toml              the shipped table merged with the best finds
    logs/workerK.log          what each worker did; NEW and IMPROVED lines for finds
```

## Merging finds

`<domain>/<table>.toml` is in the shipped format.

1. For the degrees that changed, copy the entries into `src/data`.
2. Run
   ```sh
   julia --project -t auto scripts/certify_tables.jl triangle:35,54,55 tetrahedron:23,24 cube:13
   ```
   It refines each seed to 160 bits, checks that it is the correctly rounded value, and records
   the residual.
3. `scripts/verify_tables.jl` then checks exactness, sharpness and positivity independently.
4. Update the provenance note in `src/data/PROVENANCE.toml` for the new degrees.

## Code

| File | Contents |
|---|---|
| `campaign.jl` | one worker: options, then `run_worker` |
| `engine.jl` | the pool, the five tasks, the positive-weight solve, candidates |
| `domains.jl` | one adapter per domain (triangle, tetrahedron, square, cube, sphere), over the package's orbit machinery |
| `launch.jl`, `status.jl` | start workers; summarise a campaign |

The engine is domain-independent. Each domain is an adapter in `domains.jl` of about 35 lines,
covering:
- its table file;
- its degree step: 1 on simplices, 2 on boxes, whose rules have odd degree only;
- its equation count;
- conversion between orbit structures and the table's format;
- the orbit-type counts used for balancing;
- wrappers for the package's fit, moment system, moves, orbit growth, chain, margins and
  100-digit refinement.

## The sphere

The package's own O_h-symmetric sphere rules stop at degree 17 and resume at 133: Lebedev's
published rules cover the degrees between, but only through Lebedev.jl, whose licence the
package cannot pass on. Everything from 19 to 131 was answered by the 6000-point degree-133
rule. `sphere:fill:1:19-131` grows in-house rules through that gap from the degree-17 rule.

Two things had kept the package's own search (`grow_and_eliminate_octahedral`) at degree 17,
where `notes/v0.4-spheres.md` records a stall at degree 19 with a residual near 3e-6:

- **The fit was ill-conditioned.** `fit_octahedral` solves the moment equations in the
  invariants p₄ᵃp₆ᵇ, whose Jacobian has condition number 6e5 at the degree-17 rule, about 1e19
  one orbit later and 1e48 at degree 133. Float64 cannot solve that. In the orthonormal
  invariant harmonics (`OctahedralHarmonicSystem`, which the 100-digit refinement already uses)
  it is 23 at degree 17 and 2.6e3 at degree 133. The sphere adapter fits in those.
- **Every rule of degree 18 or more needs a d orbit**, the only kind with points off the nine
  mirror planes of O_h. The product of the planes' linear forms, xyz(x²−y²)(y²−z²)(z²−x²), has
  degree 9 and changes sign under the group; its square is an invariant of degree 18 that is
  zero, with zero gradient, at every point of every other orbit kind, and has a positive
  integral, so no rule without a d orbit can integrate it. From the degree-17 rule, every
  growth without a d orbit stalls at the same residual (2.2 in the harmonics) and every one with
  it can converge. Lebedev's 146-point degree-19 rule is a₁, a₂, a₃, three 24-point orbits and
  one d orbit. The adapter's chain always adds a d orbit first.

A first 15-minute run from degree 17 reached degree 77, matching Lebedev's counts at degrees
19, 21 and 23 (146, 170 and 194 points) and 2–16% above them higher up, where the fill mode had
moved on as soon as each degree was near square:

| Degree | 19 | 21 | 23 | 25 | 31 | 41 | 53 | 65 | 77 |
|---|---|---|---|---|---|---|---|---|---|
| Found | 146 | 170 | 194 | 266 | 362 | 638 | 1034 | 1518 | 2094 |
| Lebedev | 146 | 170 | 194 | 230 | 350 | 590 | 974 | 1454 | 2030 |

## The newer domains

The 4-cube, disk, ball, pyramid, wedge and 4-simplex tables (main, 74b73a9) start at low
degree. Their adapters are thin: the 4-cube, disk and ball are the box adapter with `D = 4` or
the package's round geometry (`refine/box_search.jl`); the pyramid and wedge go through the
package's generic orbit search (`refine/orbit_search.jl`); the 4-simplex is the simplex adapter
with `N = 5`. A first six-minute run of each, all six at once on a busy machine, reached: wedge
7–14, disk 23–29, ball 15, pyramid 10 and 4-cube 19. The 4-simplex did not get past degree 6:
grown from the degree-6 rule, the degree-7 fits converge, but to rules outside the simplex or
with a negative weight; multistart over the structures of a point count, as
`scripts/generate_simplex4_seeds.jl` does to degree 4, is the likelier way in.

## Extending it to other domains

The engine needs nothing new. A wedge or pyramid adapter does need orbit machinery the package
doesn't have yet for those domains; at present it has only their conical and product rules.

**Wedge** (triangular prism, `T × [-1, 1]`).
- *Symmetry:* S₃ × Z₂, of order 12.
- *Orbits:* a triangle orbit (centroid, vertex type or general) times a height of either
  `z = 0` or `±z`, so six orbit types in all.
- *Invariant basis:* the invariant polynomials are products of S₃-invariant polynomials on the
  triangle and even polynomials in `z`. An orthonormal invariant basis is therefore the
  triangle's `invariant_basis` times even orthonormal Legendre polynomials, with no new
  projection to compute.
- *Rest:* the moment system, margins and moves follow the box code, with one more kind of
  move: sending `z` to zero.

**Pyramid** (square base, apex on the axis).
- *Symmetry:* C₄ᵥ, of order 8, acting on each horizontal slice as on the square.
- *Orbits:* a square orbit of the slice (centre, `(a, 0)`, `(a, a)` or `(a, b)`, scaled by
  `1 − z`) at a height `z`.
- *Invariant basis:* this has to be computed as on the simplex: the Reynolds projection of an
  orthonormal pyramid basis, block by block in degree, checked against its Molien series.
  Options for the basis are the rational Bergot–Cohen–Duruflé functions or a polynomial basis
  that is exact on the pyramid's polynomial space.

**Both** also need:
- first rules from multistart at low degree, as `generate_box_seeds.jl` does;
- a refinement function (Gauss–Newton with the precision floor of `refine_seed`);
- a seed table with provenance.

This search then extends them.

## Notes

- The search state is the pool. Workers can be stopped and started at any time, and a new
  worker picks up from wherever the pool is.
- This replaces the single-process scripts `scripts/search_campaign.jl` and
  `scripts/search_box_campaign.jl`, which ran the package's chains back to back; the chain task
  runs the same chains. They are in the history up to commit 6e11464.
