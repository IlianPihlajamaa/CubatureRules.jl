# CubatureRules.jl

Quadrature and cubature rules generated on demand, at arbitrary order and arbitrary
precision.

```julia
using CubatureRules

r = rule(Simplex{2}(); degree = 20)      # construct once — visibly
integrate(x -> exp(x[1] * x[2]), r)      # apply as often as you like
```

Construction is always a separate, visible step. There is deliberately no
`integrate(f, domain; degree = 20)`: the package does not cache rules, so that overload
would quietly rebuild the rule on every call inside a loop. Build the rule once, bind it to
a variable, and reuse it.

## Precision

```julia
rule(Simplex{2}(); degree = 20, digits = 200)        # BigFloat, 200 decimal digits
rule(Simplex{2}(); degree = 20, T = BigFloat)        # the ambient BigFloat precision
rule(Simplex{2}(); degree = 9, T = Rational{BigInt}) # exact (Grundmann–Möller)
rule(Simplex{2}(); degree = 12, T = Double64)        # any AbstractFloat
```

## Choosing a rule

```julia
available(Simplex{2}(); degree = 17)      # ranked candidates, nothing constructed
compare(Simplex{2}(), 17)                 # side by side, including filtered-out ones
rule(Simplex{2}(); degree = 17, positive = true, interior = true)
rule(GrundmannMöller(), Simplex{4}(); degree = 11, T = Rational{BigInt})   # explicit family
```

## Hot loops and meshes

```julia
s = static(rule(Simplex{2}(); degree = 6))    # isbits, SVector-backed
integrate(f, s)                               # allocation-free
integrate(f, s, triangle)                     # affine map applied on the fly
integrate(f, s, cells)                        # sum over a vector of simplices
```

When `integrate` is called from your own function that takes the integrand as an argument,
annotate it as `f::F ... where {F}`. Julia does not specialise on a function argument
that is only passed through, and without specialisation the call dispatches dynamically.

## Checking a rule

```julia
check(r)               # exactness, sharpness, weight sum, interior, positivity, symmetry
certificate(r)         # how refinement went: residual, cond(J), guard digits
provenance(r)          # where it came from, and why the selector chose it
cite(r)                # BibTeX (or style = :apa, :plain)
```
