# Exactness claims

Each rule states which functions it integrates exactly, through its exactness claim, a
subtype of [`ExactnessClaim`](@ref):

| Claim | Meaning |
|---|---|
| [`PolynomialDegree`](@ref)`(d)` | exact for all polynomials of total degree `≤ d` on the domain |
| [`SpanOf`](@ref)`(basis)` | exact for every function in the span of `basis` |
| [`NoClaim`](@ref)`()` | exact for no particular space |

`exactness(r)` returns the claim, and `degree(r)` returns `d` for a `PolynomialDegree` claim.

## Why there is more than one kind

Most families have a polynomial degree. Some do not: tanh-sinh and the other
double-exponential rules converge very quickly for analytic integrands, even with endpoint
singularities, but they are not exact for any polynomial space. Giving them a degree would
be wrong, and leaving them out would lose a useful tool. `NoClaim` lets them exist in the
package without pretending.

Families with `NoClaim` are never offered in answer to a `degree` request. They are checked
by a convergence study rather than by exact integration; see [Verification](verification.md).

## Tensor products

A tensor product of one-dimensional rules is exact on a tensor-product polynomial space,
which is larger than the total-degree space. A 4 × 4 Gauss rule integrates `x⁶y⁶` exactly,
but its claim is total degree 7. The selector compares candidates by total degree, so that
is the claim reported.

## What operations do to a claim

Operations that change a rule state what they do to its claim:

| Operation | Effect on `PolynomialDegree(d)` |
|---|---|
| `map_to(r, dom)` (affine) | preserved |
| `subdivide(r, dom, n)` | preserved |
| `static(r)` | preserved |
| `transform(r, φ, Jφ)` | replaced by `NoClaim`, unless the caller passes `claim` |
| `duffy(r; power)` | replaced by `NoClaim` |

A polynomial composed with an affine map is a polynomial of the same degree, so the first
three keep the claim. A nonlinear map does not preserve polynomials, so the claim is
removed rather than reduced. If a caller knows that the claim holds for a particular map,
it can be passed to `transform`; the provenance then records that the claim was asserted by
the caller and not derived by the package.

## Sharpness

A claim of degree `d` is checked in both directions: the rule must be exact at degree `d`
and not exact at degree `d + 1`. Without the second check, a rule could report a lower
degree than it has, and the selector, which ranks rules by points for a given degree, would
make worse choices.
