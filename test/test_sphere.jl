using CubatureRules, Test, LinearAlgebra
const CR = CubatureRules

@testset "sphere domain" begin
    @test CR.dimension(Sphere{3}()) == 3            # the coordinate count, not the surface dimension
    @test isreference(Sphere{3}()) && !isreference(Sphere((1, 0, 0), 2))
    @test CR.reference(Sphere((1, 0, 0), 2)) == Sphere{3}()
    @test Float64(measure(Sphere{2}())) ≈ 2π
    @test Float64(measure(Sphere{3}())) ≈ 4π
    @test Float64(measure(Sphere{4}())) ≈ 2π^2                  # 2π^{D/2}/Γ(D/2)
    @test Float64(measure(Sphere((0, 0, 0), 3))) ≈ 4π * 9       # r^{D-1}, the surface Jacobian
    @test indomain([1.0, 0.0, 0.0], Sphere{3}())
    @test !indomain([1.0, 1.0, 0.0], Sphere{3}())
    @test indomain([0.6, 0.8, 0.0], Sphere{3}())
    # a sphere has no boundary, so every point of it is interior — the registry's filter
    # would otherwise reject every rule on one
    @test isinterior([0.0, 0.0, 1.0], Sphere{3}())
    @test indomain([1.0, 2.0, 3.0] .+ [0.0, 0.0, 2.5], Sphere((1.0, 2.0, 3.0), 2.5))
    @test sprint(show, Sphere{3}()) == "Sphere{3}()"
    @test occursin("Sphere((1", sprint(show, Sphere((1.0, 0.0, 0.0), 2.0)))
    @test_throws ArgumentError Sphere((0, 0, 0), 0)
    @test eltype(CR.convert_domain(Float64, Sphere{3}()).centre) === Float64
end

@testset "sphere moments" begin
    @test CR.sphere_moment(3, (0, 0, 0)) ≈ measure(Sphere{3}())
    @test CR.sphere_moment(3, (1, 0, 0)) == 0                   # any odd exponent kills it
    @test CR.sphere_moment(3, (2, 1, 0)) == 0
    @test Float64(CR.sphere_moment(3, (2, 0, 0))) ≈ 4π / 3
    @test Float64(CR.sphere_moment(3, (4, 0, 0))) ≈ 4π / 5
    @test Float64(CR.sphere_moment(3, (2, 2, 0))) ≈ 4π / 15
    @test Float64(CR.sphere_moment(2, (2, 0))) ≈ π
    # x² + y² + z² = 1 on the sphere, so those three moments must sum to the area
    @test sum(CR.sphere_moment(3, α) for α in ((2, 0, 0), (0, 2, 0), (0, 0, 2))) ≈ measure(Sphere{3}())
    @test_throws ArgumentError CR.sphere_moment(3, (1, 1))
end

@testset "spherical harmonics" begin
    # orthonormality, tested through a rule exact to twice the harmonic degree
    n = 4
    r = rule(Sphere{3}(); degree = 2n + 1)
    Y = zeros(BigFloat, CR.harmonic_length(n))
    P = zeros(BigFloat, CR.legendre_length(n))
    G = zeros(BigFloat, CR.harmonic_length(n), CR.harmonic_length(n))
    for i in 1:npoints(r)
        CR.real_harmonics!(Y, P, n, nodes(r)[i])
        G .+= BigFloat(weights(r)[i]) .* Y .* Y'
    end
    @test maximum(abs, G - I) < 1e-12
    # the indexing is the degree-graded one
    @test CR.harmonic_block(0) == 1:1 && CR.harmonic_block(2) == 5:9
    @test CR.harmonic_length(3) == 16
    # the circle basis is orthonormal in the same sense
    rc = rule(Sphere{2}(); degree = 9)
    F = zeros(BigFloat, CR.circle_length(4))
    Gc = zeros(BigFloat, length(F), length(F))
    for i in 1:npoints(rc)
        CR.fourier_circle!(F, 4, nodes(rc)[i])
        Gc .+= BigFloat(weights(rc)[i]) .* F .* F'
    end
    @test maximum(abs, Gc - I) < 1e-12
end

@testset "SphereProduct" begin
    for d in (0, 1, 3, 8, 12)
        r = rule(Sphere{3}(); degree = d)
        v = check(r)
        @test family(r) == "SphereProduct"
        @test npoints(r) == cld(d + 1, 2) * (d + 1)
        @test degree(r) == d
        @test v.exact && v.sharp === true && v.positive && v.interior
        @test sum(weights(r)) ≈ 4π rtol = 1e-13
        @test all(x -> indomain(x, Sphere{3}()), nodes(r))
    end
    for d in (1, 4, 9)
        r = rule(Sphere{2}(); degree = d)
        v = check(r)
        @test npoints(r) == d + 1                  # optimal on the circle
        @test v.exact && v.sharp === true && v.positive
        @test sum(weights(r)) ≈ 2π rtol = 1e-13
    end
    # integrating something with a known value
    r = rule(Sphere{3}(); degree = 6)
    @test integrate(x -> x[3]^2, r) ≈ 4π / 3 rtol = 1e-13
    @test integrate(x -> x[1]^2 * x[2]^2, r) ≈ 4π / 15 rtol = 1e-13
    @test integrate(x -> x[1] * x[2] * x[3], r) ≈ 0 atol = 1e-14
    # arbitrary precision
    rb = rule(Sphere{3}(); degree = 15, digits = 50)
    @test passed(check(rb))
    @test abs(sum(weights(rb)) - 4big(π)) < big(10.0)^-48
    @test_throws NoRuleError rule(Sphere{3}(); degree = 3, T = Rational{BigInt})
    # the registry knows about it, and says what it cannot do
    @test first(available(Sphere{3}(); degree = 5)).family == "SphereProduct"
    @test isempty(available(Sphere{5}(); degree = 5))
    msg = try rule(Sphere{5}(); degree = 5) catch e; sprint(showerror, e) end
    @test occursin("Sphere{5}", msg)
end

@testset "spheres off the reference" begin
    s = rule(Sphere((1.0, 2.0, 3.0), 2.5); degree = 9)
    @test npoints(s) == 50
    @test sum(weights(s)) ≈ 4π * 2.5^2 rtol = 1e-13
    @test all(x -> indomain(x, domain(s)), nodes(s))
    @test exactness(s) == PolynomialDegree(9)        # a similarity preserves the claim
    @test passed(check(s))
    @test occursin("similarity image", last(provenance(s).path))
    @test integrate(x -> 1.0, s) ≈ 4π * 2.5^2 rtol = 1e-12
    # the mean of |x - c|² over the sphere is r²
    @test integrate(x -> sum(abs2, x .- [1.0, 2.0, 3.0]), s) / (4π * 2.5^2) ≈ 6.25 rtol = 1e-12
    @test_throws ArgumentError map_to(rule(Sphere{3}(); degree = 3), Simplex{3}())
end

@testset "domains still to come" begin
    @test_throws CR.NotYetImplemented Ball{3,Float64}()
    @test_throws CR.NotYetImplemented Pyramid()
    # a sphere in a dimension whose harmonics are not implemented says so
    @test_throws ArgumentError CR.verification_basis(Sphere{5}(), 3, BigFloat, false)
end
