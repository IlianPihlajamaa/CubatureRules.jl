# Pooled search for fully symmetric rules

A long-running search for new and smaller seeded rules on the triangle, tetrahedron, square
and cube. Many single-threaded workers share one pool of exact rules on disk, instead of each
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
repeatedly picks one of four tasks:

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

**Choosing a task.** Apart from the 25% of chains, a worker chooses by what the pool holds at
its target degree:
- *grow* while the pool is nearly empty;
- *prune* while nothing there is near square;
- otherwise mostly *finish*, with some grow and prune to keep the pool varied.

**Candidates.** Any degree-`n` rule with fewer points than the best known, at any excess, is a
candidate:
- it is refined to 100 digits by the package (`refine_symmetric`, `refine_box`);
- it is checked to be positive and interior;
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
| `engine.jl` | the pool, the four tasks, the positive-weight solve, candidates |
| `domains.jl` | one adapter per domain, over the package's orbit machinery |
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

## Extending it to wedges and pyramids

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
- a refinement function (Gauss–Newton with the precision floor of `refine_box`);
- a seed table with provenance.

This search then extends them.

## Notes

- The search state is the pool. Workers can be stopped and started at any time, and a new
  worker picks up from wherever the pool is.
- The older single-process scripts, `scripts/search_campaign.jl` and
  `scripts/search_box_campaign.jl`, remain. They run the package's chains directly.
