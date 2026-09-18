using CubatureRules, Test, Aqua

# GenericLinearAlgebra adds generic methods to LinearAlgebra functions (a known, accepted
# piracy of that dependency), so piracy checks are limited to this package's own methods.
Aqua.test_all(CubatureRules; ambiguities = false, piracies = (; treat_as_own = []), persistent_tasks = false)
@test isempty(Test.detect_ambiguities(CubatureRules; recursive = true))
