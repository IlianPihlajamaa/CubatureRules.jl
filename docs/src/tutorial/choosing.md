# I want to choose a rule

When you call `rule(domain; degree)`, the package chooses a family for you. To see which
families it considered, use [`available`](@ref). It lists the candidates without building
any of them:

```@repl choosing
using CubatureRules
available(Simplex{2}(); degree = 17)
```

The list is in the order the selector ranks them: fewest points first, then families that
are computed from a formula before families that start from a stored table, then by name.
`rule` builds the first candidate that meets your requirements, and stores the ranking in
the rule's provenance:

```@repl choosing
provenance(rule(Simplex{2}(); degree = 17)).selection
```

[`compare`](@ref) gives the same list with two extra columns: whether the family can produce
exact rational rules, and a rough estimate of how expensive it is to construct.

```@repl choosing
compare(Simplex{2}(), 17)
```

## Positive weights and interior nodes

```@repl choosing
r = rule(Simplex{2}(); degree = 17, positive = true, interior = true);
family(r), npoints(r)
```

With `positive = true`, families that can produce negative weights are left out. With
`interior = true`, families that can place nodes on the boundary are left out. If no family
remains, `rule` throws a [`NoRuleError`](@ref) explaining what was excluded; see
[Troubleshooting](troubleshooting.md).

Negative weights can cost accuracy through cancellation, and can give a negative result for
a positive integrand, so it is usually a good idea to ask for positive weights.

## Choosing the family yourself

Pass the family as the first argument:

```@repl choosing
rule(GrundmannMöller(), Simplex{4}(); degree = 11)
```

The selector is skipped in that case, and the provenance notes that you chose the family:

```@repl choosing
provenance(rule(ConicalProduct(), Simplex{2}(); degree = 5)).selection
```

Some families are built out of other families. `ConicalProduct(GaussJacobi())` builds a
simplex rule from one-dimensional Gauss–Jacobi rules, `TensorProduct(Fejer(), 2)` builds a
rule on a square from Fejér rules, and `BallProduct(SphereProduct())` combines a radial rule
with a rule on the sphere. When you let the selector choose, it considers all of these
combinations along with the other families.

## Asking for a number of points

Some families can be asked for a number of points instead of a degree:

```@repl choosing
r = rule(GaussLegendre(), Interval(); npoints = 12);
degree(r)
```

This only works when you name the family, because a number of points only has a meaning for
a specific construction.

## About the degree

`rule` returns a rule of *at least* the degree you asked for. Many families only exist at
some degrees, so you may get a higher one:

```@repl choosing
degree(rule(Simplex{2}(); degree = 8))
```

A few families have no polynomial degree at all. `TanhSinh`, `ExpSinh` and `SinhSinh`
converge quickly for many integrands, including some with endpoint singularities, but they
are not exact on any space of polynomials. Their exactness is recorded as
[`NoClaim`](@ref), and `rule(domain; degree)` never selects them. Ask for them by name, as in
`rule(TanhSinh(5), Interval())`. See [Exactness claims](../design/claims.md) for how the
package deals with this.

## Next steps

- [I want arbitrary precision](precision.md)
- [I want to verify a rule](verifying.md)
- [Selection](../design/selection.md) describes the ranking in full, and how families are
  found.
