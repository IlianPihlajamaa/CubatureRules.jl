# API

The rule families are documented on the [Rule families](families.md) page.

```@autodocs
Modules = [CubatureRules]
Private = false
Filter = t -> !(t in (XiaoGimbutas, FullySymmetric, GrundmannMöller, ConicalProduct, TensorProduct, NewtonCotes, Fejer, TanhSinh, GaussJacobi, GaussLegendre, GaussLaguerre, GaussHermite, ExpSinh, SinhSinh, SphereProduct, Lebedev, BallProduct))
```
