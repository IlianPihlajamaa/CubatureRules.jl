# Rule families

| Family | Domain | Degrees | Derivation | Points | Output types |
|---|---|---|---|---|---|
| `XiaoGimbutas` | triangle | 1–26 | seeded | minimal (Xiao–Gimbutas counts) | floating |
| `FullySymmetric` | tetrahedron | 1–9 | seeded | smallest found by in-house search | floating |
| `GrundmannMöller` | any simplex | any (odd) | derived | ``\sum_{i=0}^{s} \binom{s-i+D}{D}`` | floating, `Rational{BigInt}` |
| `ConicalProduct` | any simplex | any | derived | ``\lceil (d+1)/2 \rceil^D`` | floating |
| `GaussJacobi` | interval, Jacobi weight | any | derived | ``\lceil (d+1)/2 \rceil`` | floating |

## XiaoGimbutas

```@docs
XiaoGimbutas
```

Fully symmetric, positive and interior triangle rules at the minimal point counts of
H. Xiao and Z. Gimbutas, *Comput. Math. Appl.* 59 (2010) 663–676,
doi:10.1016/j.camwa.2009.10.027.

The seeds are MIT-licensed. They were generated in-house from orbit structures alone, and
no published numbers were used. Minimal rules are not unique; the shipped rule is the one
with the largest minimum barycentric coordinate among those found.

## FullySymmetric

```@docs
FullySymmetric
```

Fully symmetric (S₄), positive and interior tetrahedron rules. Unlike the triangle table,
nothing here comes from a paper. The point count at each degree is the smallest that the
in-house search reached, walking upwards from a count below which no fully symmetric rule
can exist. The counts are not proven minimal, and the family is not named after any
published one.

## GrundmannMöller

```@docs
GrundmannMöller
```

A. Grundmann and H. M. Möller, *SIAM J. Numer. Anal.* 15 (1978) 282–290,
doi:10.1137/0715019.

## ConicalProduct

```@docs
ConicalProduct
```

A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall, 1971.

## GaussJacobi

```@docs
GaussJacobi
GaussLegendre
```

Seeds come from Golub–Welsch (G. H. Golub and J. H. Welsch, *Math. Comp.* 23 (1969)
221–230). They are refined by Newton on the three-term recurrence, and the weights come
from the Christoffel function.
