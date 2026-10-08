using CubatureRules, LinearAlgebra, Test
const CR = CubatureRules

# Fully symmetric rules on the disk, the ball, the pyramid, the wedge and the 4-simplex: the
# moment systems against rules known in closed form, their Jacobians, a search at low degree,
# and the shipped tables.

"The Jacobian of `sys` at `θ` against central differences."
function jacobian_error(sys, θ)
    _, J = sys(θ)
    Jfd = hcat([(sys(θ + 1e-6 * (1:length(θ) .== k); jacobian = false)[1] -
                 sys(θ - 1e-6 * (1:length(θ) .== k); jacobian = false)[1]) / 2e-6 for k in eachindex(θ)]...)
    return maximum(abs, J - Jfd)
end

"Every shipped entry of at most `maxpts` points refines to a rule that verifies, with the claimed symmetry."
function check_shipped(fam, dom, entries; maxpts = 120)
    for e in filter(e -> e.npoints <= maxpts, entries)
        r = rule(fam, dom; degree = e.degree)
        @test npoints(r) == e.npoints && degree(r) == e.degree
        v = check(r)
        @test passed(v) && v.positive && v.interior && v.symmetric === true
    end
end

@testset "disk and ball" begin
    # four points at radius √½ (degree 3 on the disk), six at √(3/5) (degree 3 on the ball)
    r, _ = CR.RoundMomentSystem(CR.BoxStructure([[1]], 2), 3, Float64)([π / 4, sqrt(0.5)])
    @test maximum(abs, r) < 1e-14
    r, _ = CR.RoundMomentSystem(CR.BoxStructure([[1]], 3), 3, Float64)([2π / 9, sqrt(0.6)])
    @test maximum(abs, r) < 1e-14
    for D in (2, 3)
        s = CR.BoxStructure(D == 2 ? [[], [1], [2], [1, 1]] : [[], [1], [2], [1, 1, 1]], D)
        θ = D == 2 ? [0.3, 0.2, 0.5, 0.1, 0.4, 0.05, 0.6, 0.3] : [0.3, 0.2, 0.5, 0.1, 0.4, 0.05, 0.6, 0.3, 0.1]
        @test jacobian_error(CR.RoundMomentSystem(s, 7, Float64), θ) < 1e-7
    end
    @test length(CR.round_entries(2)) >= 6 && length(CR.round_entries(3)) >= 4
    check_shipped(FullySymmetric(), Disk(), CR.round_entries(2))
    check_shipped(FullySymmetric(), Ball{3}(), CR.round_entries(3))
    # chosen over the product of a radial and an angular rule
    @test family(rule(Disk(); degree = 7)) == "FullySymmetric"
    @test family(rule(Ball{3}(); degree = 5)) == "FullySymmetric"
    @test passed(verify(rule(FullySymmetric(), Disk(); degree = 7, digits = 40)))
end

@testset "pyramid" begin
    # the basis is orthonormal, and only its constant has a nonzero integral
    n = 4
    pr = rule(ConicalProduct(), Pyramid(); degree = 2n)
    terms = CR.pyramid_terms(n)
    basis_at(x) = begin
        z = x[3]
        pξ, _ = CR.legendre_orthonormal(n, x[1] / (1 - z))
        pη, _ = CR.legendre_orthonormal(n, x[2] / (1 - z))
        map(terms) do (i, j, k)
            a = 2(i + j) + 2
            q, _ = CR.jacobi_orthonormal_d(k, a, 0, 2z - 1)
            S = i == j ? pξ[i + 1] * pη[i + 1] : (pξ[i + 1] * pη[j + 1] + pξ[j + 1] * pη[i + 1]) / sqrt(2)
            S * (1 - z)^(i + j) * sqrt(2.0^(a + 1)) * q[k + 1]
        end
    end
    G = sum(w * basis_at(x) * transpose(basis_at(x)) for (x, w) in zip(nodes(pr), weights(pr)))
    @test maximum(abs, G - I) < 1e-13
    s = CR.PyramidStructure([[], [1], [2], [1, 1]])
    θ = [0.1, 0.3, 0.2, 0.4, 0.15, 0.6, 0.05, 0.5, 0.3, 0.02, 0.7, 0.2]
    @test jacobian_error(CR.PyramidMomentSystem(s, 7, Float64), θ) < 1e-7
    # a search at low degree finds Witherden & Vincent's point counts
    @test [CR.npoints(CR.orbit_multistart(CR.PyramidStructure, d; npts = 1:12, nstarts = 16)[1]) for d in 1:3] == [1, 5, 6]
    @test !isempty(CR.pyramid_entries())
    check_shipped(FullySymmetric(), Pyramid(), CR.pyramid_entries(); maxpts = 50)
    @test family(rule(Pyramid(); degree = 5)) == "FullySymmetric"
    @test passed(verify(rule(FullySymmetric(), Pyramid(); degree = 4, digits = 40)))
    @test CR.monomial_moment(Pyramid(), (0, 0, 0)) == 4 // 3
    @test CR.monomial_moment(Pyramid(), (2, 0, 1)) == 2 // 45                           # ∫ x² z
end

@testset "wedge" begin
    # the centroid at z = 0 is the degree-1 rule
    r, _ = CR.WedgeMomentSystem(CR.WedgeStructure([[3, 0]]), 1, Float64)([1.0])
    @test maximum(abs, r) < 1e-14
    s = CR.WedgeStructure([[3, 0], [3, 2], [2, 1, 0], [2, 1, 2], [1, 1, 1, 0], [1, 1, 1, 2]])
    θ = [0.05, 0.03, 0.4, 0.02, 0.2, 0.01, 0.15, 0.6, 0.02, 0.1, 0.3, 0.015, 0.25, 0.1, 0.7]
    @test jacobian_error(CR.WedgeMomentSystem(s, 6, Float64), θ) < 1e-7
    @test CR.n_wedge_equations(6) == CR.n_equations(CR.WedgeMomentSystem(s, 6, Float64))
    @test [CR.npoints(CR.orbit_multistart(CR.WedgeStructure, d; npts = 1:12, nstarts = 16, slack = 4)[1])
           for d in 1:3] == [1, 5, 8]
    @test !isempty(CR.wedge_entries())
    check_shipped(FullySymmetric(), Wedge(), CR.wedge_entries(); maxpts = 50)
    @test family(rule(Wedge(); degree = 5)) == "FullySymmetric"
    @test passed(verify(rule(FullySymmetric(), Wedge(); degree = 4, digits = 40)))
    @test CR.monomial_moment(Wedge(), (1, 0, 2)) == big(1) // 6 * 2 // 3
end

@testset "4-simplex" begin
    n = 3
    X, W, _, _ = CR.conical_work(4, n + 1, 128)
    b = CR.SimplexBasis{4,Float64}(n)
    G = zeros(CR.basis_length(b), CR.basis_length(b))
    for (x, w) in zip(X, W)
        φ, _ = CR.evaluate!(b, Float64.(x); gradient = false)
        G .+= Float64(w) .* φ .* transpose(φ)
    end
    @test maximum(abs, G - I) < 1e-13
    # the general kernel is the tetrahedral basis at D = 3
    bt = CR.SimplexBasis{3,Float64}(4)
    φg, Gg = zeros(CR.basis_length(bt)), zeros(CR.basis_length(bt), 3)
    CR.simplex_dubiner!(φg, Gg, CR.SimplexWorkspace{Float64}(3, 4), [0.2, 0.3, 0.1])
    φt, Gt = CR.evaluate!(bt, [0.2, 0.3, 0.1])
    @test maximum(abs, φg - φt) < 1e-12 && maximum(abs, Gg - Gt) < 1e-11
    # the invariant basis has the dimensions of the Molien series
    @test CR.invariant_basis(5, 4).ranks == [CR.molien_coefficient(5, k) for k in 0:4] == [1, 0, 1, 1, 2]
    s = CR.SymmetricStructure([[5], [4, 1], [3, 2], [2, 2, 1]], 5)
    θ = [0.01, 0.005, 0.1, 0.004, 0.15, 0.002, 0.1, 0.25]
    @test jacobian_error(CR.SymmetricMomentSystem(s, 4, Float64), θ) < 1e-7
    @test !isempty(CR.simplex_entries(4))
    check_shipped(FullySymmetric(), Simplex{4}(), CR.simplex_entries(4); maxpts = 60)
    # fewer points than the other positive rule; Grundmann–Möller, with 7 points and a negative
    # weight, stays the default at degree 3
    @test npoints(rule(FullySymmetric(), Simplex{4}(); degree = 3)) < npoints(rule(ConicalProduct(), Simplex{4}(); degree = 3))
end

@testset "symmetry checks of the new groups" begin
    for (dom, d) in ((Pyramid(), 4), (Wedge(), 4))
        r = rule(FullySymmetric(), dom; degree = d)
        g = CR.properties(FullySymmetric(), dom, d).symmetry
        xs, ws = nodes(r), weights(r)
        @test CR.check_symmetry(g, xs, ws, 1e-12) === true
        ys = copy(xs)
        ys[end] = ys[end] + CR.SVector(1e-6, 0.0, 0.0)
        @test CR.check_symmetry(g, ys, ws, 1e-12) === false
    end
end
