using CubatureRules, Test
import CubatureRules: isreference, candidates   # public, not exported
const CR = CubatureRules

@testset "R^D and the Gaussian space" begin
    @test RealLine() === RealSpace{1}()          # the line is the one-dimensional case
    @test CR.dimension(RealSpace{4}()) == 4
    @test isreference(RealSpace{4}()) && measure(RealSpace{4}()) == Inf
    @test indomain([1.0, -2.0, 3.0], RealSpace{3}())
    @test !indomain([1.0, Inf], RealSpace{2}())
    @test sprint(show, RealSpace{1}()) == "RealLine()"
    @test sprint(show, RealSpace{4}()) == "RealSpace{4}()"
    # the weight is radial, so it reads the same in any dimension
    @test CR.GaussianWeight()(2.0) ≈ exp(-4)
    @test CR.GaussianWeight()([1.0, 1.0]) ≈ exp(-2)
    for D in 1:4
        @test Float64(measure(GaussianSpace(D))) ≈ π^(D / 2)
    end
    @test GaussianSpace(1) == HermiteLine()
end

@testset "Gaussian moments" begin
    @test Float64(CR.gaussian_moment(3, (0, 0, 0))) ≈ π^1.5
    @test CR.gaussian_moment(3, (1, 0, 0)) == 0         # any odd exponent kills it
    @test CR.gaussian_moment(2, (2, 1)) == 0
    @test Float64(CR.gaussian_moment(3, (2, 0, 0))) ≈ π^1.5 / 2
    @test Float64(CR.gaussian_moment(3, (4, 0, 0))) ≈ 3π^1.5 / 4
    @test Float64(CR.gaussian_moment(3, (2, 2, 0))) ≈ π^1.5 / 4
    # a product of one-dimensional moments, which is where it comes from
    @test CR.gaussian_moment(2, (2, 4)) ≈ CR.gaussian_moment(1, (2,)) * CR.gaussian_moment(1, (4,))
    @test_throws ArgumentError CR.gaussian_moment(3, (2, 2))
end

@testset "spheres in any dimension" begin
    # the areas of S^1 … S^4, written out rather than recomputed from the formula under test
    areas = Dict(2 => 2π, 3 => 4π, 4 => 2π^2, 5 => 8π^2 / 3)
    for D in 2:5
        @test Float64(measure(Sphere{D}())) ≈ areas[D]
        for d in (1, 3, 6)
            r = rule(SphereProduct(), Sphere{D}(); degree = d)
            v = check(r)
            n, m = CR.sphere_product_shape(Sphere{D}(), d)
            @test npoints(r) == n^(D - 2) * m
            @test v.exact && v.sharp === true && v.positive
            @test sum(weights(r)) ≈ Float64(measure(Sphere{D}())) rtol = 1e-12
            @test all(x -> indomain(x, Sphere{D}()), nodes(r))
        end
    end
    # the circle stays optimal, and S² keeps its old shape
    @test npoints(SphereProduct(), Sphere{2}(), 7) == 8
    @test npoints(SphereProduct(), Sphere{3}(), 7) == 4 * 8
    # balls come along for free, since BallProduct takes any sphere family
    r = rule(BallProduct(SphereProduct()), Ball{4}(); degree = 5)
    @test passed(check(r))
    @test sum(weights(r)) ≈ Float64(measure(Ball{4}())) rtol = 1e-12
end

@testset "GaussianProduct: Stroud's E_n^{r²}, derived" begin
    for D in 2:4
        g = GaussianSpace(D)
        for d in (1, 3, 5, 9)
            r = rule(GaussianProduct(), g; degree = d)
            v = check(r)
            @test degree(r) >= d
            @test v.exact && v.positive
            @test sum(weights(r)) ≈ Float64(measure(g)) rtol = 1e-12
            @test npoints(r) == CR.gaussian_radial_points(d) * npoints(SphereProduct(), Sphere{D}(), d)
        end
    end
    # the radial factor is cheap because odd degrees never reach it
    @test CR.gaussian_radial_points(1) == 1
    @test CR.gaussian_radial_points(9) == 3
    # known integrals in three dimensions
    r = rule(GaussianProduct(), GaussianSpace(3); degree = 6)
    @test integrate(x -> x[1]^2, r) ≈ π^1.5 / 2 rtol = 1e-12
    @test integrate(x -> x[1]^2 * x[2]^2, r) ≈ π^1.5 / 4 rtol = 1e-12
    @test integrate(x -> x[1]^4, r) ≈ 3π^1.5 / 4 rtol = 1e-12
    @test integrate(x -> x[1] * x[2], r) ≈ 0 atol = 1e-14      # odd, so it vanishes
    # arbitrary precision
    rb = rule(GaussianProduct(), GaussianSpace(3); degree = 7, digits = 40)
    @test passed(check(rb))
    dev = CR.with_bits(() -> abs(sum(weights(rb)) - measure(GaussianSpace(3))), CR.digits_to_bits(40) + 64)
    @test dev < big(10.0)^-37
    @test_throws NoRuleError rule(GaussianProduct(), GaussianSpace(3); degree = 3, T = Rational{BigInt})
    # the one-dimensional case belongs to Gauss–Hermite, which is cheaper
    @test isempty(candidates(GaussianProduct, GaussianSpace(1), PolynomialDegree(5)))
    @test family(rule(HermiteLine(); degree = 5)) == "GaussHermite"
    # and the selector finds it in higher dimensions, pairing it with the cheapest angular
    # family it can — which on S² is Lebedev, not the product rule
    names = [row.family for row in available(GaussianSpace(3); degree = 5)]
    @test "GaussianProduct(SphereProduct)" in names
    @test first(names) == "GaussianProduct(Lebedev)"
    rl = rule(GaussianSpace(3); degree = 5)
    @test passed(check(rl))
    @test npoints(rl) < npoints(rule(GaussianProduct(SphereProduct()), GaussianSpace(3); degree = 5))
end
