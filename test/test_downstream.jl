# A family defined in a separately precompiled package is found by the selector with no
# registration step (PLAN §2.5). This is the guard against `subtypes` ever being cached in
# a `const`: that would be evaluated when CubatureRules was precompiled, before
# TestFamilyPackage existed.
using CubatureRules, Test

# test/downstream is an implicit environment (a directory of packages): its packages resolve
# their dependencies through the rest of the load path, i.e. the test environment.
push!(LOAD_PATH, joinpath(@__DIR__, "downstream"))
try
    @eval using TestFamilyPackage
finally
    pop!(LOAD_PATH)
end
const CentroidRule = TestFamilyPackage.CentroidRule

@test CentroidRule in families()
@test candidates(CentroidRule, Simplex{2}(), PolynomialDegree(1)) == [CentroidRule()]
a = available(Simplex{2}(); degree = 1)
@test first(a).family == "CentroidRule"
r = rule(Simplex{2}(); degree = 1)
@test family(r) == "CentroidRule"
@test npoints(r) == 1
@test passed(check(r))
@test occursin("CentroidRule (1 points)", provenance(r).selection)
@test family(rule(Simplex{3}(); degree = 0)) == "CentroidRule"
# degree 2 is out of its range, so selection falls back to the package's families
@test family(rule(Simplex{2}(); degree = 2)) != "CentroidRule"
