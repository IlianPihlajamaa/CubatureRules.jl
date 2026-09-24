using CubatureRules
using Test

const CR = CubatureRules

@testset "CubatureRules.jl" begin
    @testset "core types and domains" include("test_core.jl")
    @testset "symmetry and orbit algebra" include("test_symmetry.jl")
    @testset "refinement" include("test_refine.jl")
    @testset "families" include("test_families.jl")
    @testset "delegated families" include("test_delegated.jl")
    @testset "moment-defined weights" include("test_moments.jl")
    @testset "unbounded domains" include("test_unbounded.jl")
    @testset "spheres" include("test_sphere.jl")
    @testset "octahedral orbit algebra" include("test_octahedral.jl")
    @testset "octahedral grow and eliminate" include("test_octahedral_grow.jl")
    @testset "Lebedev" include("test_lebedev.jl")
    @testset "Lebedev.jl interop" include("test_upstream_lebedev.jl")
    @testset "balls" include("test_ball.jl")
    @testset "Gaussian space and spheres in any dimension" include("test_gaussian.jl")
    @testset "registry and selection" include("test_registry.jl")
    @testset "application" include("test_apply.jl")
    @testset "claim preservation" include("test_transport.jl")
    @testset "composition" include("test_composition.jl")
    @testset "embedded rules" include("test_embedded.jl")
    @testset "verification" include("test_verify.jl")
    @testset "number types" include("test_numbertypes.jl")
    @testset "presentation" include("test_presentation.jl")
    @testset "provenance manifest" include("test_provenance.jl")
    @testset "reproducibility" include("test_reproducibility.jl")
    @testset "property-based" include("test_properties.jl")
    @testset "downstream family discovery" include("test_downstream.jl")
    @testset "scripts" include("test_scripts.jl")
    @testset "Aqua" include("test_aqua.jl")
end
