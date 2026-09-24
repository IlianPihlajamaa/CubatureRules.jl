# Interval

```@setup interval
using CubatureRules, Markdown
function family_table(dom)
    fmt(r) = last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"
    rows = ["| Family | Degrees | Derivation |", "|---|---|---|"]
    for row in available(dom)
        der = row.derivation isa CubatureRules.Derived ? "derived" : "seeded"
        push!(rows, "| `$(row.family)` | $(fmt(row.degrees)) | $der |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

## Domains

```@docs
Interval
WeightedDomain
JacobiWeight
```

On `Interval()`:

```@example interval
family_table(Interval())
```

## Gauss–Jacobi

```@docs
GaussJacobi
GaussLegendre
```

The `n`-point rule for the weight `(1-x)^α (1+x)^β` has degree `2n - 1`, the highest
possible. The seed comes from the Golub–Welsch method: the eigenvalues of the Jacobi matrix
of the three-term recurrence, in `Float64`. Newton's method on the recurrence refines the
nodes to working precision, and the weights come from the Christoffel function
`wᵢ = μ₀ / Σₖ pₖ(xᵢ)²`.

References: G. H. Golub and J. H. Welsch, *Math. Comp.* 23 (1969) 221–230; G. Szegő,
*Orthogonal Polynomials*, AMS Colloquium Publications 23 (1939).

## Newton–Cotes

```@docs
NewtonCotes
```

The weights are computed in exact rational arithmetic, so `T = Rational{BigInt}` is
supported. Newton–Cotes rules are mainly useful when equally spaced nodes are required, or
as a source of exact rules.

## Fejér

```@docs
Fejer
```

Fejér rules have positive weights at every order and do not use the endpoints. Their degree
is lower than that of a Gauss rule with the same number of points, but the nodes and weights
are given by simple formulas.

## Gauss–Patterson

```@docs
GaussPatterson
```

Each level keeps the nodes of the one below, so the pair gives an error estimate from a single
set of function values:

```@example interval
e = embedded(GaussPatterson(), Interval(); degree = 20)
integrate(x -> 1 / (1 + 25x^2), e; error = true)
```

Every level is derived here from the one below: the new nodes are the roots of a polynomial
found from one linear system, bracketed between the old nodes, and the weights are
interpolatory. The rules become ill-conditioned to compute: the 255-point rule needs about 43
extra digits of working precision and the 511-point rule about 95. The reason is that each
rule comes very close to integrating one degree more than it claims. The 127-point rule misses
degree 192 by only `1e-20`, which is below `Float64` rounding. [`verify`](@ref) reports this
as sharpness that cannot be resolved at this precision, quoting the miss measured during
construction. At 40 digits the miss is resolved:

```@example interval
verify(rule(GaussPatterson(), Interval(); npoints = 127, digits = 40))
```

In `Float64` the 255-point rule takes about 4 s to build and the 511-point rule about 20 s.

## Tanh-sinh

```@docs
TanhSinh
```

Tanh-sinh rules have no polynomial degree, so they report [`NoClaim`](@ref) and are
requested by level: `rule(TanhSinh(5), Interval())`. They are checked by a convergence study
instead of exact integration:

```@example interval
seq = [rule(TanhSinh(m), Interval()) for m in 2:6];
CubatureRules.verify_convergence(seq, x -> 1 / sqrt(1 - x^2), Float64(π))
```

They are useful for integrands with endpoint singularities. On `1/√(1-x²)` with 201 points,
tanh-sinh reaches an error of `5e-8`, where Gauss–Legendre reaches `9e-3`.

The accuracy is limited by how well `1 - x` can be represented near the endpoints. The
outermost node lies a few rounding units from the endpoint, so an integrand evaluated as a
function of `x` loses about half the working digits there: about `3e-8` in `Float64` and
`3e-26` at 50 digits. Asking for more digits, or writing the integrand in terms of the
distance to the endpoint, avoids this.

## Delegated families

These families are computed by other packages; see
[External providers](../design/providers.md).

```@docs
GaussKronrod
Lobatto
Radau
ClenshawCurtis
```

## Weights given by moments

```@docs
ModifiedChebyshev
MomentWeight
OrdinaryMoments
LogWeight
CubatureRules.christoffel
```

For a weight without a classical family, `ModifiedChebyshev` computes the Gauss rule from
the weight's moments with Wheeler's algorithm. See
[I want an unusual weight](../tutorial/weights.md) for examples and
[Measures given by moments](../design/moments.md) for the method and its conditioning.
