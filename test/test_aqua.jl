using CubatureRules, Test, Aqua

Aqua.test_all(CubatureRules; ambiguities = false, persistent_tasks = false)
# MultiFloats.jl defines `(::Type{T})(x::MultiFloat{T,N})` for every `T`, which
# detect_ambiguities pairs with every one-argument constructor here, although no
# `MultiFloat{T}` exists for those `T`; the ambiguities of this package's own methods are none
own = filter(Test.detect_ambiguities(CubatureRules; recursive = true)) do (m1, m2)
    nameof(m1.module) !== :MultiFloats && nameof(m2.module) !== :MultiFloats
end
@test isempty(own)
