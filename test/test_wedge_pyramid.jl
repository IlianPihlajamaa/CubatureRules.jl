using CubatureRules, Test
const CR = CubatureRules

# Wedges and pyramids: the domains, product and conical rules on them, and integration over
# meshes of them, checked against exact monomial integrals and against cubes split into
# wedges or pyramids.

wedge_moment(a, b, c) = factorial(big(a)) * factorial(big(b)) // factorial(big(a + b + 2)) * (iseven(c) ? 2 // (c + 1) : 0 // 1)
pyramid_moment(a, b, c) = (isodd(a) || isodd(b)) ? 0 // 1 :
    4 // ((a + 1) * (b + 1)) * factorial(big(c)) * factorial(big(a + b + 2)) // factorial(big(a + b + c + 3))
function monomial_error(r, moment, d)
    m = big(0.0)
    for a in 0:d, b in 0:(d - a), c in 0:(d - a - b)
        v = setprecision(() -> integrate(x -> big(x[1])^a * big(x[2])^b * big(x[3])^c, r), BigFloat, 300)
        m = max(m, abs(v - moment(a, b, c)))
    end
    return m
end

@testset "wedge and pyramid domains" begin
    @test measure(Wedge()) == 1 && measure(Pyramid()) == 4 // 3
    @test CR.isreference(Wedge()) && CR.isreference(Pyramid())
    @test sprint(show, Wedge()) == "Wedge()" && sprint(show, Pyramid()) == "Pyramid()"
    W = Wedge((0, 0, 0), (2, 0, 0), (0, 1, 0), (1, 1, 3), (3, 1, 3), (1, 2, 3))
    P = Pyramid((0, 0, 0), (2, 0, 0), (3, 1, 0), (1, 1, 0), (1, 1, 2))
    @test measure(W) == 3 && measure(P) == 4 // 3       # base area × height; base area × height / 3
    @test !CR.isreference(W) && occursin("Wedge((0, 0, 0)", sprint(show, W))
    @test indomain((1, 2 // 3, 1), W) && isinterior((1, 2 // 3, 1), W)       # the centroid of the slice z = 1
    @test !indomain((1 // 2, 1 // 4, 1), W)                                  # below that slice's lower edge
    @test indomain((0, 0, 0), W) && !isinterior((0, 0, 0), W) && !indomain((5, 5, 5), W)
    @test isinterior((3 // 2, 1 // 2, 1 // 2), P) && !indomain((0, 0, 1), P)
    # only affine images of the reference shapes
    @test_throws ArgumentError Wedge((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 2), (0, 1, 1))
    @test_throws ArgumentError Pyramid((0, 0, 0), (1, 0, 0), (1, 2, 0), (0, 1, 0), (0, 0, 1))
    @test_throws ArgumentError Wedge((0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 0, 1), (1, 0, 1), (2, 0, 1))
    @test_throws ArgumentError Pyramid((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 2, 0))
    @test_throws ArgumentError Wedge((0, 0, 0), (1, 0, 0))
    # a float wedge from a mesh, translated to rounding, is accepted
    @test measure(Wedge((0.1, 0.2, 0.3), (1.1, 0.2, 0.3), (0.1, 1.2, 0.3), (0.1, 0.2, 1.3), (1.1, 0.2, 1.3), (0.1, 1.2, 1.3))) ≈ 0.5
end

@testset "wedge rules" begin
    for d in 0:12
        r = rule(Wedge(); degree = d)
        @test degree(r) >= d && family(r) == "WedgeProduct"
        v = verify(r)
        @test passed(v) && v.positive && v.interior
    end
    @test monomial_error(rule(Wedge(); degree = 8), wedge_moment, 8) < 1e-14
    r = rule(Wedge(); degree = 10, digits = 40)
    @test passed(verify(r)) && monomial_error(r, wedge_moment, 10) < 1e-39
    # the triangle factor follows the triangle table, and the conical rule beyond it
    @test npoints(rule(Wedge(); degree = 10)) == npoints(rule(XiaoGimbutas(), Simplex{2}(); degree = 10)) * 6
    dmax = last(CR.degree_range(XiaoGimbutas(), Simplex{2}()))
    r = rule(Wedge(); degree = dmax + 1)
    @test occursin("ConicalProduct", provenance(r).path[2]) && degree(r) >= dmax + 1
    @test monomial_error(r, wedge_moment, 6) < 1e-13
    @test npoints(WedgeProduct(), Wedge(), 10) == npoints(rule(Wedge(); degree = 10))
end

@testset "pyramid rules" begin
    for d in 0:12
        r = rule(Pyramid(); degree = d)
        @test degree(r) >= d && family(r) == "ConicalProduct" && npoints(r) == cld(d + 1, 2)^3 || d == 0
        v = verify(r)
        @test passed(v) && v.positive && v.interior
    end
    @test monomial_error(rule(Pyramid(); degree = 8), pyramid_moment, 8) < 1e-14
    r = rule(Pyramid(); degree = 11, digits = 40)
    @test passed(verify(r)) && monomial_error(r, pyramid_moment, 11) < 1e-39
end

@testset "meshes of wedges and pyramids" begin
    f(x) = x[1]^2 * x[2] + x[3]^3 + x[1] * x[2] * x[3] + 1
    # the unit cube as two wedges over the triangles of its bottom face
    W1 = Wedge((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1))
    W2 = Wedge((1, 0, 0), (1, 1, 0), (0, 1, 0), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    exact = 1 // 6 + 1 // 4 + 1 // 8 + 1
    r = rule(Wedge(); degree = 3)
    @test integrate(f, r, [W1, W2]) ≈ exact rtol = 1e-14
    @test integrate(f, map_to(r, W1)) + integrate(f, map_to(r, W2)) ≈ exact rtol = 1e-14
    @test passed(verify(map_to(r, W1)))
    # the cube [-1, 1]³ as six pyramids, one on each face, apex at the centre
    o = (0, 0, 0)
    faces = [((-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1)), ((-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)),
             ((-1, -1, -1), (1, -1, -1), (1, -1, 1), (-1, -1, 1)), ((-1, 1, -1), (1, 1, -1), (1, 1, 1), (-1, 1, 1)),
             ((-1, -1, -1), (-1, 1, -1), (-1, 1, 1), (-1, -1, 1)), ((1, -1, -1), (1, 1, -1), (1, 1, 1), (1, -1, 1))]
    pyramids = [Pyramid(q..., o) for q in faces]
    @test sum(measure, pyramids) == 8
    g(x) = x[1]^2 * x[2]^2 + x[3]^4 + 3x[1] * x[3] + 2
    r = rule(Pyramid(); degree = 4)
    @test integrate(g, r, pyramids) ≈ 8 // 9 + 8 // 5 + 16 rtol = 1e-14
    @test passed(verify(map_to(r, pyramids[3])))
end

@testset "selection and listing" begin
    @test [c.family for c in available(Pyramid())] == ["ConicalProduct(GaussJacobi)"] ||
          any(c -> occursin("ConicalProduct", c.family), available(Pyramid()))
    @test any(c -> occursin("WedgeProduct", c.family), available(Wedge(); degree = 6))
    @test_throws NoRuleError rule(Wedge(); degree = 3, T = Rational{BigInt})
end
