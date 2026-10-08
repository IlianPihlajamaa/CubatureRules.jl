using CubatureRules, Test
const CR = CubatureRules

# Fully symmetric rules on the square and the cube: orbits of the signed permutations, the
# moment system on symmetrised Legendre products, and the shipped tables.

@testset "box orbits and the moment system" begin
    @test [o.mult for o in CR.box_orbit_types(2)] == [Int[], [1], [2], [1, 1]]
    @test length(CR.box_orbit_types(3)) == 7
    @test [CR.orbit_size(o) for o in CR.box_orbit_types(3)] == [1, 6, 8, 12, 24, 24, 48]
    s = CR.BoxStructure([[], [1], [2, 1], [1, 1, 1]], 3)
    θ = [0.3, 0.2, 0.5, 0.1, 0.7, 0.2, 0.05, 0.8, 0.5, 0.3]
    xs, ws = CR.expand(s, θ)
    @test length(xs) == CR.npoints(s) == 1 + 6 + 24 + 48
    # the 3×3×3 Gauss rule is four orbits of the cube, and satisfies the degree-5 equations
    a, w0, w1 = sqrt(3 / 5), 8 / 9, 5 / 9
    g = CR.BoxStructure([[], [1], [2], [3]], 3)
    r, _ = CR.BoxMomentSystem(g, 5, Float64)([w0^3, w0^2 * w1, a, w0 * w1^2, a, w1^3, a])
    @test maximum(abs, r) < 1e-14
    r, _ = CR.BoxMomentSystem(g, 7, Float64)([w0^3, w0^2 * w1, a, w0 * w1^2, a, w1^3, a])
    @test maximum(abs, r) > 1e-3                                  # and not the degree-7 ones
    # the Jacobian against central differences
    sys = CR.BoxMomentSystem(s, 7, Float64)
    _, J = sys(θ)
    Jfd = hcat([(sys(θ + 1e-6 * (1:length(θ) .== k))[1] - sys(θ - 1e-6 * (1:length(θ) .== k))[1]) / 2e-6
                for k in eachindex(θ)]...)
    @test maximum(abs, J - Jfd) < 1e-8
    @test CR.n_equations(CR.BoxMomentSystem(s, 6, Float64)) == CR.n_equations(CR.BoxMomentSystem(s, 7, Float64))
end

# Checking a rule costs about its point count squared: 1.6 s for the 120-point square rule of
# degree 25, 47 s for the 1280-point cube rule of degree 29 and 7 minutes for the 3548-point one
# of degree 41. Every entry is checked in the full sweep (CUBATURERULES_FULL_SWEEP=1); per
# commit, as for the simplex tables (test_families.jl), every one of at most 120 points, every
# fourth of the rest and the largest, among the entries of at most 650 points.
function box_sweep(es)
    get(ENV, "CUBATURERULES_FULL_SWEEP", "0") == "1" && return es
    affordable = filter(e -> e.npoints <= 650, sort(es; by = e -> e.degree))
    large = filter(e -> e.npoints > 120, affordable)
    keep = Set(e.degree for e in vcat(filter(e -> e.npoints <= 120, affordable), large[1:4:end], last(affordable)))
    return filter(e -> e.degree in keep, affordable)
end

@testset "shipped square, cube and 4-cube rules" begin
    for D in (2, 3, 4)
        es = CR.box_entries(D)
        @test length(es) >= (D == 4 ? 4 : 6)
        for e in box_sweep(es)
            r = rule(FullySymmetric(), Orthotope{D}(); degree = e.degree)
            @test npoints(r) == e.npoints && degree(r) == e.degree
            v = check(r)
            @test passed(v) && v.positive && v.interior && v.symmetric === true
        end
    end
    # Float64 is served from the certified table as stored
    @test certificate(rule(Orthotope{2}(); degree = 9)).iterations == 0
    # more digits than the table holds: refined, and still verified
    r = rule(FullySymmetric(), Orthotope{2}(); degree = 9, digits = 40)
    @test passed(verify(r)) && occursin("Gauss–Newton", join(provenance(r).path, " "))
    @test passed(verify(rule(FullySymmetric(), Orthotope{3}(); degree = 7, digits = 40)))
    # an even degree is served by the odd rule above it
    @test degree(rule(FullySymmetric(), Orthotope{2}(); degree = 8)) == 9
    # chosen over the tensor rule on total degree
    r = rule(Orthotope{3}(); degree = 7)
    @test family(r) == "FullySymmetric" && npoints(r) == 34
    # a mapped box, against the exact integral
    box = Orthotope((0.0, 1.0, -1.0), (2.0, 2.0, 3.0))
    r = rule(box; degree = 6)
    @test passed(check(r))
    @test integrate(x -> x[1]^2 * x[2]^3 + x[3]^4, r) ≈ 8 / 3 * 15 / 4 * 4 + 2 * 1 * 244 / 5 rtol = 1e-13
end
