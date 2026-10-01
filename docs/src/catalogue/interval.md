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

## Weights given by a function

```@docs
StieltjesDiscretization
FunctionWeight
PointMass
```

```@example interval
w = FunctionWeight(x -> exp(-x), 0, 2) + PointMass(1, 1 // 2);
r = rule(WeightedDomain(Interval(0, 2), w); degree = 15, digits = 40);
last(provenance(r).path, 4)
```

Written as function weights, the Legendre, Jacobi, Laguerre and Hermite weights reproduce
the classical rules to the last digit, and the weight above gives the same 40 digits as
`ModifiedChebyshev` on its exact moments. Rules on a function weight are verified against
the weight's moments, computed by a separate discretisation at the verification precision.

The discrete Stieltjes procedure is implemented here. PolyChaos.jl has one, but it computes
in `Float64` whatever its input. QuadGK builds Gauss rules for a weight function at any
precision, but evaluates Chebyshev series at every node on every step and was 30 times
slower on these sizes.

## Oscillatory weights

```@docs
Filon
Oscillatory
```

For `∫ f(x) e^{iωx} dx` with large `ω`, a rule for the weight `e^{iωx}` needs only enough
points for `f`, not for the oscillation. Here `∫₀¹ eˣ e^{1000ix} dx` with 11 points:

```@example interval
r = rule(WeightedDomain(Interval(0, 1), Oscillatory(1000)); degree = 10, digits = 30);
z = 1 + 1000im;
npoints(r), integrate(exp, r), (exp(big(z)) - 1) / z
```

The rules are interpolatory on the Gauss–Lobatto points. With the end points among the
nodes, the error at a fixed number of points falls like `ω⁻²` (A. Iserles and S. P. Nørsett,
*Proc. R. Soc. A* 461 (2005) 1383–1399, doi:10.1098/rspa.2004.1401):

```@example interval
for ω in (10, 100, 1000, 10^4)
    r = rule(WeightedDomain(Interval(0, 1), Oscillatory(ω)); degree = 8, digits = 30)
    z = 1 + ω * im
    println("ω = ", lpad(ω, 5), ":  error ", abs(integrate(exp, r) - (exp(big(z)) - 1) / z))
end
```

The weights come from the Legendre moments of the weight, `∫ P_k(x) e^{iκx} dx = 2iᵏ j_k(κ)`
on `[-1, 1]`, with `j_k` the spherical Bessel functions. Interpolation on the Lobatto points
needs no linear solve: the Lobatto rule integrates `P_j P_k` exactly except for `j = k = N`,
so the interpolant's Legendre coefficients are discrete inner products. Nothing cancels as
`ω → 0`, unlike the classical closed-form Filon weights, which are differences of large
terms there. The `j_k` are computed upwards from the closed forms of `j₀` and `j₁`, at a
precision raised by the digits that recurrence loses above `k ≈ κ`. Verification computes
them downwards by Miller's algorithm, so the two share only `sin` and `cos`.

For `e^{iωx}` the weights are complex and the nodes real; `Oscillatory(ω, cos)` and
`Oscillatory(ω, sin)` give real rules. On an interval centred at 0, `cos` and `sin` have a
parity that the symmetric rule shares, and the rule is then exact one degree higher. As
`ω → 0` the rule tends to the Gauss–Lobatto rule, and misses degree `N + 1` by an amount
too small for a low-precision rule to show; the certificate records the miss, as for
Gauss–Patterson.

## Principal values and finite parts

```@docs
SingularGauss
PrincipalValue
FinitePart
```

The Cauchy principal value `⨍ f(x) w(x) / (x − t) dx` and the Hadamard finite part
`⨎ f(x) w(x) / (x − t)² dx`, with a Jacobi weight `w`, are linear functionals of `f`, and
a rule can be exact for them on polynomials like any other. Here the finite part of
`√(1 − x²) / (x − 1/3)²`, which is `−π`:

```@example interval
dom = WeightedDomain(Interval(-1, 1), FinitePart(1 // 3; α = 1 // 2, β = 1 // 2));
r = rule(dom; degree = 20, digits = 30);
npoints(r), degree(r), integrate(x -> one(x), r)
```

```@example interval
last(provenance(r).path, 3)
```

Everything is computed from the Hilbert transform of the weight, `ρ₀(t) = ⨍ w / (x − t)`,
and its derivative `σ₀(t)`. For a Jacobi weight these have a closed form in a
hypergeometric series (W. Gautschi and J. Wimp, *BIT* 27 (1987) 203–215), with a
logarithmic series for integer exponents. They are evaluated at two precisions, with
more guard bits until the two agree, because an exponent close to an integer makes two
large terms cancel. The integrals `ρ_k(t)` and `σ_k(t)` of the orthogonal polynomials
follow from the three-term recurrence.

The principal-value rule is D. B. Hunter's, *Numer. Math.* 19 (1972) 419–424,
doi:10.1007/BF01404924: the `n` Gauss–Jacobi nodes and `t`, exact to degree `2n`. The
finite-part rule also has `n` nodes and `t` and is exact to degree `2n`. Writing a
polynomial as `a + b(x − t) + (x − t)² p(x)`, the last term is integrated exactly when the
`n` nodes are those of a rule for `w` exact to degree `2n − 2`. With `q_k` the orthogonal
polynomials of `w`, those nodes are the zeros of `q_n − c q_{n−1}`, and requiring the term
`b(x − t)` to come out right fixes `c = ρ_n(t)/ρ_{n−1}(t)`. The published Gauss-type
finite-part rules either need `f'(t)` or are exact only to degree `n − 1`. These zeros are
real but can leave the interval, which happens for about half of all `n`: when
`t = cos θ`, the bad values of `n` come in runs of length about `π / (2θ)`. The
construction then tries larger `n`, up to twice the minimum `n₀`. If none of those works,
it uses the interpolatory rule on `2n₀` Gauss nodes and `t`. That rule is exact to the
requested degree and has no more points than the largest finite-part rule tried.

The weights are signed, and `Σ|wᵢ|` grows as `t` approaches a node, so the construction
takes one more point when `t` lies within an eighth of a gap of a node. A rule that comes
out exact beyond its construction's degree, as a rule centred in a symmetric weight
does on every odd polynomial, claims the higher degree.

These rules are verified against the functional applied to the orthogonal polynomials of
the interval. Those integrals are computed by a Gauss–Jacobi rule applied to divided
differences, without the recurrence the construction uses, so the two share only `ρ₀` and
`σ₀`. The tests check those two against classical closed forms, and the rules' integrals
of `eˣ` against tanh-sinh at 90 digits on the subtracted integrand. At 40 digits, both
rules agree to within `10⁻³⁸`, for non-integer, integer and nearly integer exponents, and
for `t` within `10⁻²⁰` of an end.
