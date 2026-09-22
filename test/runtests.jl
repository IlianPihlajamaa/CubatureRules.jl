using CubatureRules
using Test

const CR = CubatureRules

@testset "CubatureRules.jl" begin
    @testset "core types and domains" include("test_core.jl")
    @testset "symmetry and orbit algebra" include("test_symmetry.jl")
    @testset "refinement" include("test_refine.jl")
    @testset "families" include("test_families.jl")
    @testset "delegated families" include("test_delegated.jl")
    @testset "unbounded domains" include("test_unbounded.jl")
    @testset "spheres" include("test_sphere.jl")
    @testset "registry and selection" include("test_registry.jl")
    @testset "application" include("test_apply.jl")
    @testset "claim preservation" include("test_transport.jl")
    @testset "composition" include("test_composition.jl")
    @testset "sequences and tolerance" include("test_sequence.jl")
    @testset "verification" include("test_verify.jl")
    @testset "number types" include("test_numbertypes.jl")
    @testset "presentation" include("test_presentation.jl")
    @testset "provenance manifest" include("test_provenance.jl")
    @testset "reproducibility" include("test_reproducibility.jl")
    @testset "property-based" include("test_properties.jl")
    @testset "downstream family discovery" include("test_downstream.jl")
    @testset "Aqua" include("test_aqua.jl")
end
