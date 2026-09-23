# Families whose numbers come from upstream packages (PLAN §0.1). QuadGK is a hard
# dependency; QuadratureRules.jl is weak, so its families only work once it is loaded —
# which is also a second, independent check that `subtypes` discovery works across a
# package boundary.
using CubatureRules, Test
const CR = CubatureRules

@testset "Gauss–Kronrod (QuadGK)" begin
    for d in (2, 5, 11, 20)
        r = rule(GaussKronrod(), Interval(); degree = d)
        v = check(r)
        @test isodd(npoints(r))
        @test degree(r) >= d
        @test sum(weights(r)) ≈ 2
        @test v.exact && v.sharp === true
        @test v.positive && v.interior && v.symmetric === true
    end
    # the degree is 3n+1, or 3n+2 for odd n
    @test CR.kronrod_degree(2) == 7 && CR.kronrod_degree(3) == 11
    @test npoints(GaussKronrod(), Interval(), 7) == 2 * CR.kronrod_halves(7) + 1
    r = rule(GaussKronrod(), Interval(); degree = 20, digits = 60)
    @test passed(check(r))
    @test occursin("QuadGK", provenance(r).seed_source)
    @test provenance(r).citations[1].key == "Johnson_QuadGK"
    @test_throws NoRuleError rule(GaussKronrod(), Interval(); degree = 5, T = Rational{BigInt})
end

# Loading it is all that is needed; the extension does the rest.
#
# `import`, not `using`: QuadratureRules exports `nodes`, `weights` and `QuadratureRule`
# too, and every test file is included into `Main`, so a `using` here would make those
# names ambiguous for every file that follows. Julia 1.11 turned that into an error at
# use rather than a warning at import.
import QuadratureRules

@testset "QuadratureRules.jl families, once it is loaded" begin
    for f in (Lobatto(), Radau(:right), Radau(:left), ClenshawCurtis())
        @test CR.missing_dependency(f) === nothing        # the extension is loaded in the test env
        for d in (2, 5, 11)
            r = rule(f, Interval(); degree = d)
            v = check(r)
            @test degree(r) >= d
            @test sum(weights(r)) ≈ 2
            @test v.exact && v.sharp === true && v.positive
            @test occursin("QuadratureRules", provenance(r).seed_source)
            @test provenance(r).citations[1].key == "Kraus_QuadratureRules"
        end
    end
    # endpoints are nodes, so these are not interior rules, and the selector can filter them
    @test !check(rule(Lobatto(), Interval(); degree = 5)).interior
    @test !check(rule(ClenshawCurtis(), Interval(); degree = 5)).interior
    a = available(Interval(); degree = 7)
    names = [row.family for row in a]
    @test "Lobatto" in names && "ClenshawCurtis" in names && "Radau(:right)" in names
    @test first(names) == "GaussLegendre"                  # still the fewest points
    ai = available(Interval(); degree = 7, interior = true)
    @test !("Lobatto" in [row.family for row in ai])
    # arbitrary precision comes through the delegation unchanged
    r = rule(Lobatto(), Interval(); degree = 11, digits = 50)
    @test precision(first(weights(r))) == CR.digits_to_bits(50)
    @test passed(check(r))
    @test abs(sum(weights(r)) - 2) < big(10.0)^-48
    # Lobatto endpoints are exactly ±1
    @test CR.nodes(r)[1] == -1 && CR.nodes(r)[end] == 1
end
