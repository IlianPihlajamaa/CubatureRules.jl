using CubatureRules, Test
const CR = CubatureRules

@testset "ball domain" begin
    @test Float64(measure(Disk())) ≈ π
    @test Float64(measure(Ball{3}())) ≈ 4π / 3
    @test Float64(measure(Ball{4}())) ≈ π^2 / 2
    @test Float64(measure(Ball((0, 0, 0), 2))) ≈ 4π / 3 * 8            # r^D, the volume Jacobian
    @test Disk() == Ball{2}()
    @test isreference(Ball{3}()) && !isreference(Ball((1, 0, 0), 2))
    @test CR.reference(Ball((1, 0, 0), 2)) == Ball{3}()
    # a ball has an inside, unlike a sphere
    @test indomain([0.0, 0.0, 0.0], Ball{3}()) && isinterior([0.0, 0.0, 0.0], Ball{3}())
    @test indomain([1.0, 0.0, 0.0], Ball{3}())                          # the boundary belongs to it
    @test !isinterior([1.0, 0.0, 0.0], Ball{3}())                       # but is not interior
    @test !indomain([1.0, 1.0, 0.0], Ball{3}())
    @test isinterior([1.5, 2.0, 3.0], Ball((1.0, 2.0, 3.0), 2.5))
    @test sprint(show, Ball{3}()) == "Ball{3}()"
    @test occursin("Ball((1", sprint(show, Ball((1.0, 0.0, 0.0), 2.0)))
    @test_throws ArgumentError Ball((0, 0, 0), 0)
end

@testset "ball moments" begin
    @test CR.ball_moment(3, (0, 0, 0)) ≈ measure(Ball{3}())
    @test CR.ball_moment(3, (1, 0, 0)) == 0
    @test CR.ball_moment(3, (2, 1, 0)) == 0
    @test Float64(CR.ball_moment(3, (2, 0, 0))) ≈ 4π / 15               # sphere's 4π/3 over 5
    @test Float64(CR.ball_moment(2, (2, 0))) ≈ π / 4
    # the defining relation: the sphere's moment over |α| + D
    for α in ((0, 0, 0), (2, 0, 0), (2, 2, 0), (4, 2, 2))
        @test CR.ball_moment(3, α) ≈ CR.sphere_moment(3, α) / (sum(α) + 3)
    end
end

@testset "BallProduct" begin
    for d in (0, 1, 3, 6, 9)
        r = rule(BallProduct(SphereProduct()), Ball{3}(); degree = d)
        v = check(r)
        @test degree(r) >= d
        @test npoints(r) == CR.ball_radial_points(d) * CR.npoints(SphereProduct(), Sphere{3}(), d)
        @test v.exact && v.positive && v.interior && v.weights_sum_ok
        @test sum(weights(r)) ≈ 4π / 3 rtol = 1e-13
        @test all(x -> indomain(x, Ball{3}()), nodes(r))
    end
    # the angular factor can be any sphere family, and Lebedev is the cheap one
    r = rule(BallProduct(Lebedev()), Ball{3}(); degree = 9)
    @test passed(check(r))
    @test npoints(r) == 5 * 38
    @test occursin("Lebedev", CR.describe_family(BallProduct(Lebedev())))
    @test npoints(r) < npoints(rule(BallProduct(SphereProduct()), Ball{3}(); degree = 9))
    # the selector ranks the combinations for you
    @test first(available(Ball{3}(); degree = 9)).family == "BallProduct(Lebedev)"
    # the disk
    rd = rule(Disk(); degree = 5)
    @test passed(check(rd))
    @test sum(weights(rd)) ≈ π rtol = 1e-13
    @test integrate(x -> x[1]^2, rd) ≈ π / 4 rtol = 1e-12
    # known integrals on the 3-ball
    r = rule(Ball{3}(); degree = 6)
    @test integrate(x -> x[3]^2, r) ≈ 4π / 15 rtol = 1e-12
    @test integrate(x -> sum(abs2, x), r) ≈ 4π / 5 rtol = 1e-12
    @test integrate(x -> x[1] * x[2], r) ≈ 0 atol = 1e-14
    # arbitrary precision
    rb = rule(Ball{3}(); degree = 7, digits = 40)
    @test passed(check(rb))
    dev = CR.with_bits(() -> abs(sum(weights(rb)) - measure(Ball{3}())), CR.digits_to_bits(40) + 64)
    @test dev < big(10.0)^-38
    @test_throws NoRuleError rule(Ball{3}(); degree = 3, T = Rational{BigInt})
end

@testset "balls off the reference" begin
    b = Ball((1.0, 2.0, 3.0), 2.5)
    s = rule(BallProduct(Lebedev()), b; degree = 7)
    @test sum(weights(s)) ≈ 4π / 3 * 2.5^3 rtol = 1e-12
    @test all(x -> indomain(x, b), nodes(s))
    @test exactness(s) == PolynomialDegree(7)
    @test passed(check(s))
    @test occursin("similarity image", last(provenance(s).path))
    @test integrate(x -> 1.0, s) ≈ 4π / 3 * 2.5^3 rtol = 1e-12
    # the mean of |x - c|² over a ball of radius R is 3R²/5
    @test integrate(x -> sum(abs2, x .- [1.0, 2.0, 3.0]), s) / (4π / 3 * 2.5^3) ≈ 3 * 2.5^2 / 5 rtol = 1e-12
end
