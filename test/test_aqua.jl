using CubatureRules, Test, Aqua

Aqua.test_all(CubatureRules; ambiguities = false, persistent_tasks = false)
@test isempty(Test.detect_ambiguities(CubatureRules; recursive = true))
