# I want arbitrary precision

Pass `digits` to get a rule in `BigFloat` with that many decimal digits:

```@repl prec
using CubatureRules
r = rule(Simplex{2}(); degree = 20, digits = 200);
nodes(r)[1][1]
```

The nodes and weights are computed at this precision. They are not a `Float64` table
converted to `BigFloat`.

## Choosing the number type

```@repl prec
rule(Interval(); degree = 9, digits = 50);             # BigFloat, 50 decimal digits
rule(Interval(); degree = 9, T = BigFloat);            # the current BigFloat precision
rule(Interval(); degree = 9, T = Float32);             # any AbstractFloat
rule(Simplex{2}(); degree = 9, T = Rational{BigInt});  # exact, where possible
```

- `digits = d` gives a `BigFloat` rule with `d` decimal digits.
- `T = BigFloat` uses the current BigFloat precision, so it follows an enclosing
  `setprecision`.
- `T` can be any `AbstractFloat`, including `Double64` from DoubleFloats.jl or `Float128`
  from Quadmath.jl. The rule is computed in `BigFloat` and converted to `T` at the end.
- `T = Rational{BigInt}` asks for an exact rule.

## Exact rules

Some families have rational nodes and weights, and can return them exactly:

```@repl prec
r = rule(Simplex{2}(); degree = 5, T = Rational{BigInt});
family(r), nodes(r)[1]
```

If no family can produce the type you ask for, `rule` throws a [`NoRuleError`](@ref) rather
than rounding. On an interval, for example, the Gauss–Legendre nodes are irrational, so an
exact request gets a `NewtonCotes` rule instead.

## How accurate the result is

During the computation the package works with more precision than you asked for. The number
of extra digits is chosen from the condition number of the equations that define the rule,
so a harder problem gets more. The rule's [`certificate`](@ref) records the target digits,
the extra (guard) digits and the condition number:

```@repl prec
c = certificate(rule(Simplex{2}(); degree = 20, digits = 50));
c.digits, c.guard_digits, c.cond
```

All the work is done in `BigFloat`, and the result is rounded to the output type once, at the
end. Since MPFR rounds correctly, the same call gives bitwise identical results on every
platform. [Precision and guard digits](../design/precision.md) explains this in more detail.

At `Float64` and below, the minimal symmetric rules on triangles, tetrahedra and spheres
are not recomputed. The package's tables already hold them, correctly rounded and checked,
so they are returned or rounded directly. That is why their certificate shows no guard
digits and no iterations.

## Cost

At 15 digits the rule below comes straight from the table. Beyond that it is computed, which
takes longer, but not by as much as you might expect: the refinement uses Newton's method,
which roughly doubles the number of correct digits in each step.

```@repl prec
using Printf
for d in (15, 50, 200)
    t = @elapsed rule(Simplex{2}(); degree = 20, digits = d)
    @printf("%4d digits: %6.3f s\n", d, t)
end
```

Using the rule costs the same at any precision as ordinary arithmetic in that type. The extra
cost is paid once, when the rule is built.

Some constructions take minutes, such as large symmetric rules on the sphere. Pass
`verbose = true` to see their progress; see [Troubleshooting](troubleshooting.md).

## Next steps

- [I want to verify a rule](verifying.md) to check the result.
- [I want an unusual weight](weights.md) shows a case where high precision is needed to get
  an answer at all.
