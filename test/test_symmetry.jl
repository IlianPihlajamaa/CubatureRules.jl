using CubatureRules, Test
const CR = CubatureRules

@testset "Molien series" begin
    # S₃ on the triangle: 1/((1-t²)(1-t³))
    series3 = [CR.molien_coefficient(3, k) for k in 0:12]
    @test series3 == [1, 0, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3]
    # S₄ on the tetrahedron: 1/((1-t²)(1-t³)(1-t⁴))
    @test [CR.molien_coefficient(4, k) for k in 0:8] == [1, 0, 1, 1, 2, 1, 3, 2, 4]
    @test CR.n_invariants(3, 20) == 44
end

@testset "invariant basis dimension matches Molien at every degree" begin
    for n in (1, 5, 10, 20, 25)
        B = CR.invariant_basis(3, n)
        @test B.ranks == [CR.molien_coefficient(3, k) for k in 0:n]
        @test size(B.Q, 2) == CR.n_invariants(3, n)
        # orthonormal columns
        @test maximum(abs, B.Q' * B.Q - CR.LinearAlgebra.I) < 1e-12
    end
    # Q is computed with plain loops: bitwise identical on repeat
    @test CR.invariant_basis(3, 12).Q == CR.invariant_basis(3, 12).Q
end

@testset "orbit patterns" begin
    @test CR.orbit_size(CR.OrbitPattern([3], 3)) == 1
    @test CR.orbit_size(CR.OrbitPattern([2, 1], 3)) == 3
    @test CR.orbit_size(CR.OrbitPattern([1, 1, 1], 3)) == 6
    # tetrahedral patterns (v0.2 needs data, not code)
    pats = CR.orbit_pattern_types(4)
    @test [p.mult for p in pats] == [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    @test [CR.orbit_size(p) for p in pats] == [1, 4, 6, 12, 24]
    @test_throws ArgumentError CR.OrbitPattern([2, 2], 3)
end

@testset "orbit expansion is exact under the group action" begin
    s = CR.SymmetricStructure([[1, 1, 1], [2, 1], [3]], 3)
    @test CR.npoints(s) == 10
    @test CR.nunknowns(s) == 3 + 2 + 1
    θ = [0.01, 0.1, 0.25, 0.02, 0.2, 0.05]
    λs, ws = CR.expand(s, θ)
    @test length(λs) == 10
    @test all(λ -> sum(λ) ≈ 1, λs)
    for λ in λs, σ in CR.permutations_of(3)
        @test any(μ -> μ == λ[σ], λs)      # exact equality, not approximate
    end
    @test sum(ws) ≈ 6 * 0.01 + 3 * 0.02 + 0.05
end

@testset "canonical form" begin
    s = CR.SymmetricStructure([[1, 1, 1], [1, 1, 1]], 3)
    θ = [0.1, 0.3, 0.2, 0.05, 0.1, 0.15]
    θswap = [0.05, 0.15, 0.1, 0.1, 0.2, 0.3]   # same orbits, permuted and reordered
    @test CR.canonicalize(s, θ)[2] ≈ CR.canonicalize(s, θswap)[2]
end

@testset "moment system Jacobian against finite differences" begin
    e = CR.xg_entry_for(8)
    sys = CR.TriangleMomentSystem(e.structure, e.degree, Float64)
    θ = e.seed .+ 1e-3 .* sin.(1:length(e.seed))
    r, J = sys(θ)
    h = 1e-7
    Jfd = similar(J)
    for j in eachindex(θ)
        δ = zeros(length(θ)); δ[j] = h
        Jfd[:, j] = (sys(θ + δ)[1] - sys(θ - δ)[1]) / 2h
    end
    @test maximum(abs, J - Jfd) < 1e-6
    @test size(J) == (CR.n_invariants(3, e.degree), CR.nunknowns(e.structure))
end
