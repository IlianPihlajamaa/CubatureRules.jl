# Troubleshooting

## `rule` throws a `NoRuleError`

A [`NoRuleError`](@ref) means that no loaded family can produce what you asked for. The
message says why, and what would be possible instead.

Forgetting the degree gives a list of what is available on the domain:

```@repl trouble
using CubatureRules
rule(Simplex{2}())
```

Asking a specific family for a degree it does not have gives the nearest one it does have:

```@repl trouble
rule(XiaoGimbutas(), Simplex{2}(); degree = 60)
```

Requirements that no family can meet together are reported with the families that were
excluded and the reason for each:

```@repl trouble
rule(Simplex{2}(); degree = 9, positive = true, T = Rational{BigInt})
```

Here the only family with exact rational rules at this degree (`GrundmannMöller`) has
negative weights. Drop one of the two requirements, or accept a floating-point rule.

## A family needs another package

Some families use code from other packages, which have to be loaded first:

```@repl trouble
rule(Lobatto(), Interval(); degree = 9)
```

After `using QuadratureRules`, `Lobatto`, `Radau` and `ClenshawCurtis` are available and
appear in [`available`](@ref). See [External providers](../design/providers.md).

## A warning about a licence

```julia
julia> using Lebedev

julia> rule(Sphere{3}(); degree = 29)
┌ Warning: SphereProduct (450 points) was selected on Sphere{3}() at degree 29;
│ Lebedev (Lebedev.jl, GPL-3) has 302 points but its licence is not one this package
│ can pass on. Pass `copyleft = true` to use it, or
│ `CubatureRules.license_warnings!(false)` to silence this.
```

This appears when a smaller rule exists but its licence (here GPL-3, from Lebedev.jl) means
it is not chosen automatically. You get the larger rule. To use the smaller one, pass
`copyleft = true` or name the family. To turn the warning off, call
`CubatureRules.license_warnings!(false)`. Each case warns only once.

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

## Verification fails

If `check(r)` fails for a rule you built with `rule`, that is a bug in the package; please
report it with the call that produced the rule. The failing part of the report (exactness,
sharpness or structure) and the size of the residual help to narrow it down.

For a rule you built yourself, see [Verification](../design/verification.md) for what each
check tests.
