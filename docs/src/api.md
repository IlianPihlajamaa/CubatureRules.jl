# API

The rule families are documented on the [Rule families](families.md) page, and the machinery
for measures given by their moments on the [Weights given by moments](moments.md) page.

```@autodocs
Modules = [CubatureRules]
Private = false
Filter = t -> !(t in (XiaoGimbutas, FullySymmetric, GrundmannMöller, ConicalProduct, TensorProduct, NewtonCotes, Fejer, TanhSinh, GaussJacobi, GaussLegendre, GaussLaguerre, GaussHermite, ExpSinh, SinhSinh, SphereProduct, LebedevRule, UpstreamLebedev, BallProduct, ModifiedChebyshev, MomentWeight, OrdinaryMoments, CubatureRules.MonicRecurrence, CubatureRules.monic, CubatureRules.shift, CubatureRules.monomial_recurrence, CubatureRules.wheeler, CubatureRules.MomentBreakdownError))
```

Names that are `public` but not exported appear above too. Reach them as
`CubatureRules.name`, or bring them in explicitly with `using CubatureRules: name`.
