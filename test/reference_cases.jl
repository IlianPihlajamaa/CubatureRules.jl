# The rules whose hashes are pinned by the reproducibility check.
using CubatureRules

const REFERENCE_CASES = [
    "xg_triangle_d20_200digits" => () -> rule(Simplex{2}(); degree = 20, digits = 200),
    "xg_triangle_d10_float64" => () -> rule(Simplex{2}(); degree = 10),
    "xg_triangle_d15_50digits" => () -> rule(XiaoGimbutas(), Simplex{2}(); degree = 15, digits = 50),
    "gm_tet_d7_rational" => () -> rule(GrundmannMöller(), Simplex{3}(); degree = 7, T = Rational{BigInt}),
    "conical_triangle_d11_80digits" => () -> rule(ConicalProduct(), Simplex{2}(); degree = 11, digits = 80),
    "gauss_legendre_n40_100digits" => () -> rule(GaussLegendre(), Interval(); npoints = 40, digits = 100),
    "gauss_jacobi_21_float64" => () -> rule(WeightedDomain(Interval(), JacobiWeight(2, 1)); degree = 15),
]
