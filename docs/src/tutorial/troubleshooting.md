# Troubleshooting

## `rule` throws a `NoRuleError`

A [`NoRuleError`](@ref) means that no loaded family can produce what you asked for. The
message says why, and what would be possible instead.

Forgetting the degree gives a list of what is available on the domain:

```@repl trouble
using CubatureRules
rule(Simplex{2}())
```

Asking a specific family for a degree it does not have, here one beyond its table, gives the
nearest one it does have:

```@repl trouble
d = last(CubatureRules.degree_range(XiaoGimbutas(), Simplex{2}())) + 1
rule(XiaoGimbutas(), Simplex{2}(); degree = d)
```

Requirements that no family can meet together are reported with the families that were
excluded and the reason for each:

```jldoctest trouble; setup = :(using CubatureRules)
julia> rule(Simplex{2}(); degree = 9, positive = true, T = Rational{BigInt})
ERROR: No rule with positive weights and exact Rational{BigInt} nodes and weights exists on
Simplex{2}() at degree 9 among the loaded families.
  XiaoGimbutas (19 points) cannot produce Rational{BigInt}.
  ConicalProduct(GaussJacobi) (25 points) cannot produce Rational{BigInt}.
  GrundmannMöller (35 points) has negative weights.
Nearest: GrundmannMöller degree 1 (1 point, always available).
Stacktrace:
[...]
```

Here the only family with exact rational rules at this degree (`GrundmannMöller`) has
negative weights. Drop one of the two requirements, or accept a floating-point rule.

## A family needs another package

Some families use code from other packages, which have to be loaded first:

```jldoctest trouble
julia> rule(Lobatto(), Interval(); degree = 9)
ERROR: Lobatto builds its rules with QuadratureRules.jl, which is not loaded; run
`using QuadratureRules` first
Stacktrace:
[...]
```

After `using QuadratureRules`, `Lobatto`, `Radau` and `ClenshawCurtis` are available and
appear in [`available`](@ref). See [External providers](../design/providers.md).

## A warning that a smaller rule exists

At some degrees Lebedev.jl has a smaller rule on the sphere than this package's own table.
Asking for one of them without Lebedev.jl loaded says so:

```@setup trouble
# the lowest degree where Lebedev.jl has a smaller positive rule than the in-house one (its
# degree-13 rule has a negative weight, so it is never suggested)
ds = first(d for (d, n) in CubatureRules.UPSTREAM_LEBEDEV_POINTS
           if d != 13 && n < npoints(rule(LebedevRule(), Sphere{3}(); degree = d)))
CubatureRules.selection_warnings!(true)
```

```@repl trouble
rule(Sphere{3}(); degree = ds)
```

With Lebedev.jl loaded, the same request warns instead that its GPL-3 rule was passed over
for its licence, and that `copyleft = true` would allow it. These warnings appear when a
smaller rule exists that `rule` did not choose: because the package that has it is not
loaded, or because its licence (here GPL-3, from Lebedev.jl) means it is not chosen
automatically. You get the larger rule. To use the smaller one, load the
package and pass `copyleft = true`, or name the family. A rule with negative weights is
never suggested in place of a positive one. To turn the warnings off, call
`CubatureRules.selection_warnings!(false)`. Each case warns only once.

## A rule takes a long time

Most rules are built in well under a second. Large symmetric rules at high degree and high
precision can take minutes. Pass `verbose = true` to see what is happening:

```@repl trouble
rule(XiaoGimbutas(), Simplex{2}(); degree = 15, digits = 40, verbose = true);
```

The first message is printed before the expensive part starts and shows the size of the
problem. After that there is one line per iteration, with the residual, the step size, the
condition number and the time taken. `verbose = 2` shows more detail, such as the steps of
the line search.

The messages are sent with `@info`, so they can be redirected or filtered with Julia's
logging functions like any other log messages.

To stop a long construction from another task, pass a cancellation token:

```julia
tok = CubatureRules.CancellationToken()
task = Threads.@spawn rule(UpstreamLebedev(), Sphere{3}(); degree = 125, digits = 100,
                           cancel = tok)
# later
CubatureRules.cancel!(tok)          # the task stops with a CancelledError
```

On Julia 1.11, keep other `BigFloat` work out of the way while such a task runs: the
precision it sets is global there, not per task (see [Threads](../design/precision.md#Threads)).

## Verification fails

If `check(r)` fails for a rule you built with `rule`, that is a bug in the package; please
report it with the call that produced the rule. The failing part of the report (exactness,
sharpness or structure) and the size of the residual help to narrow it down.

For a rule you built yourself, see [Verification](../design/verification.md) for what each
check tests.
