using CubatureRules, Test, StaticArrays
const CR = CubatureRules

@testset "integrate is exact on polynomials up to the degree" begin
    r = rule(Simplex{2}(); degree = 12)
    for (a, b) in ((0, 0), (3, 4), (7, 5), (12, 0))
        @test integrate(x -> x[1]^a * x[2]^b, r) ≈ Float64(monomial_moment(Simplex{2}(), (a, b))) rtol = 1e-13
    end
    l = rule(Interval(); degree = 11)
    @test integrate(x -> x^10, l) ≈ 2 / 11
end

@testset "generic in the integrand's return type" begin
    r = rule(Simplex{2}(); degree = 6)
    z = integrate(x -> complex(x[1], x[2]), r)
    @test z isa ComplexF64
    @test z ≈ complex(1 / 6, 1 / 6)
    v = integrate(x -> SVector(1.0, x[1], x[2]^2), r)
    @test v isa SVector{3,Float64}
    @test v ≈ SVector(0.5, 1 / 6, 1 / 12)
    m = integrate(x -> x * x', r)
    @test m isa SMatrix{2,2,Float64}
    @test m[1, 2] ≈ 1 / 24
    # BigFloat rule, BigFloat result
    rb = rule(Simplex{2}(); degree = 8, digits = 50)
    @test integrate(x -> x[1] * x[2], rb) isa BigFloat
    @test abs(integrate(x -> x[1] * x[2], rb) - big(1) / 24) < big(10.0)^-49
end

@testset "integrate over a physical domain and a mesh" begin
    r = rule(Simplex{2}(); degree = 5)
    t = Simplex((1.0, 1.0), (3.0, 1.0), (1.0, 2.0))
    @test integrate(x -> 1.0, r, t) ≈ 1.0
    @test integrate(x -> x[1], r, t) ≈ 1.0 * (1 + 3 + 1) / 3
    @test integrate(x -> x[1]^2 * x[2], r, t) ≈ integrate(x -> x[1]^2 * x[2], map_to(r, t))
    # mesh: the unit square split into two triangles
    mesh = [Simplex((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)), Simplex((0.0, 0.0), (1.0, 1.0), (0.0, 1.0))]
    @test integrate(x -> x[1]^2 * x[2]^3, r, mesh) ≈ 1 / 12
    li = rule(Interval(); degree = 5)
    @test integrate(x -> x^2, li, Interval(0.0, 3.0)) ≈ 9.0
end

@testset "batched evaluation" begin
    r = rule(Simplex{2}(); degree = 6)
    @test integrate(X -> X[1, :] .* X[2, :], r; batch = true) ≈ integrate(x -> x[1] * x[2], r)
    l = rule(Interval(); degree = 7)
    @test integrate(X -> X .^ 2, l; batch = true) ≈ 2 / 3
    @test_throws DimensionMismatch integrate(X -> [1.0], r; batch = true)
end

# measured behind a function barrier, as in a caller's inner loop
hot_f(x) = x[1] * x[2] + sin(x[1])
hot_g(x) = exp(x)
# (`f::F` forces specialisation on the integrand, as a real kernel would)
allocs(f::F, r) where {F} = (integrate(f, r); @allocated integrate(f, r))
allocs(f::F, r, d) where {F} = (integrate(f, r, d); @allocated integrate(f, r, d))

@testset "hot path on static rules allocates nothing" begin
    s = static(rule(Simplex{2}(); degree = 6))
    @test allocs(hot_f, s) == 0
    t = Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 1.0))
    @test allocs(hot_f, s, t) == 0
    sl = static(rule(Interval(); degree = 9))
    @test allocs(hot_g, sl) == 0
    @test integrate(hot_g, sl) ≈ exp(1) - exp(-1)
    @test allocs(hot_g, sl, Interval(0.0, 2.0)) == 0
end

@testset "threaded mesh integration" begin
    r = static(rule(Simplex{2}(); degree = 5))
    m = 40
    h = 1 / m
    cells = [Simplex((i * h, j * h), ((i + 1) * h, j * h), ((i + 1) * h, (j + 1) * h)) for i in 0:(m - 1), j in 0:(m - 1)]
    cells = vcat(vec(cells),
                 vec([Simplex((i * h, j * h), ((i + 1) * h, (j + 1) * h), (i * h, (j + 1) * h)) for i in 0:(m - 1), j in 0:(m - 1)]))
    f(p) = exp(-p[1] * p[2]) * (1 + p[1]^2)
    serial = integrate(f, r, cells)
    @test integrate(f, r, cells; threaded = true) ≈ serial rtol = 1e-14
    @test serial ≈ integrate(f, rule(Simplex{2}(); degree = 12), cells) rtol = 1e-6
    # the threaded sum is over per-thread partials in thread order, so it is repeatable
    @test integrate(f, r, cells; threaded = true) == integrate(f, r, cells; threaded = true)
    @test_throws ArgumentError integrate(f, r, Simplex{2,Float64,3}[])
end
