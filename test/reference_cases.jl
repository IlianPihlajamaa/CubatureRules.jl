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
    "tensor_box3_d7_60digits" => () -> rule(TensorProduct(GaussLegendre(), 3), Orthotope{3}(); degree = 7, digits = 60),
    "symmetric_cube_d9_60digits" => () -> rule(FullySymmetric(), Orthotope{3}(); degree = 9, digits = 60),
    "symmetric_square_d15_float64" => () -> rule(FullySymmetric(), Orthotope{2}(); degree = 15),
    "newton_cotes_closed_d9_rational" => () -> rule(NewtonCotes(), Interval(); degree = 9, T = Rational{BigInt}),
    "fejer2_d15_80digits" => () -> rule(Fejer(2), Interval(); degree = 15, digits = 80),
    "kronrod_d20_60digits" => () -> rule(GaussKronrod(), Interval(); degree = 20, digits = 60),
    # the orbit families on their other domains, refined and shipped
    "symmetric_tet_d8_60digits" => () -> rule(FullySymmetric(), Simplex{3}(); degree = 8, digits = 60),
    "symmetric_simplex4_d5_60digits" => () -> rule(FullySymmetric(), Simplex{4}(); degree = 5, digits = 60),
    "symmetric_4cube_d7_float64" => () -> rule(FullySymmetric(), Orthotope{4}(); degree = 7),
    "symmetric_disk_d11_60digits" => () -> rule(FullySymmetric(), Disk(); degree = 11, digits = 60),
    "symmetric_ball_d9_60digits" => () -> rule(FullySymmetric(), Ball{3}(); degree = 9, digits = 60),
    "symmetric_wedge_d6_60digits" => () -> rule(FullySymmetric(), Wedge(); degree = 6, digits = 60),
    "symmetric_pyramid_d7_60digits" => () -> rule(FullySymmetric(), Pyramid(); degree = 7, digits = 60),
    "lebedev_d17_60digits" => () -> rule(LebedevRule(), Sphere{3}(); degree = 17, digits = 60),
    "lebedev_d151_float64" => () -> rule(LebedevRule(), Sphere{3}(); degree = 151),
    # the boundary-element kernel rules
    "duffygauss_d9_40digits" => () -> rule(WeightedDomain(SurfaceTriangle((0, 0, 0), (1, 0, 1), (3 // 10, 8 // 10, 0)),
                                                          InverseDistance((13 // 30, 8 // 30, 1 // 3))); degree = 9, digits = 40),
    "duffysinh_d5_30digits" => () -> rule(WeightedDomain(SurfaceTriangle((0, 0, 0), (1, 0, 0), (3 // 10, 8 // 10, 0)),
                                                         InverseDistance((1 // 3, 1 // 4, 1 // 100))); degree = 5, digits = 30),
    "duffyfinitepart_d4_30digits" => () -> rule(WeightedDomain(Simplex((0, 0), (1, 0), (3 // 10, 8 // 10)),
                                                               InverseDistanceGradient((2 // 5, 3 // 10), (3 // 5, 4 // 5)));
                                                degree = 4, digits = 30),
    "duffycurved_d3_float64" => () -> rule(WeightedDomain(Simplex{2}(),
                                                          CurvedInverseDistance(QuadraticTriangle((1, 0, 0), (0, 1, 0), (0, 0, 1),
                                                                                                  (1 // 2, 1 // 2, 0), (0, 1 // 2, 1 // 2),
                                                                                                  (1 // 2, 0, 1 // 2)),
                                                                                (1 // 5, 3 // 10))); degree = 3),
]
