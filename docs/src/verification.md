# Verification

```julia
v = check(rule)          # or verify(rule; degree, bits)
passed(v)
```

For a `PolynomialDegree(d)` claim, verification checks the following:

1. **Exactness** against an *orthonormal* basis up to degree `d`: Dubiner on triangles,
   Legendre or Jacobi on intervals, and barycentric monomials with exact Dirichlet moments
   for simplices of dimension 3 and up. It runs at twice the rule's precision. Each basis
   function's residual is held to the error that rounding the delivered nodes and weights
   can cause. Exact rational rules are held to zero.
2. **Sharpness**: the rule is *not* exact at degree `d + 1`, so a claimed degree cannot
   silently drift upward.
3. **Structure**: the weights sum to the measure, the nodes are interior, the weights are
   positive, and any claimed symmetry group is preserved exactly.

A rule with `NoClaim` has nothing to verify by exact integration. Use
`CubatureRules.verify_convergence` on a sequence of such rules; the result is flagged
`empirical`.

In your own test suites:

```julia
using Test
@test_exact rule 10
```
