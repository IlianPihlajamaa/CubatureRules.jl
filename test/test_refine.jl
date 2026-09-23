using CubatureRules, Test, LinearAlgebra
import CubatureRules: CancellationToken, CancelledError, ExplicitSeed, LowerDegreeSeed, cancel!   # public, not exported
const CR = CubatureRules

@testset "Gauss–Newton in generic arithmetic" begin
    # intersection of a circle and a line, solved to 300 bits
    F(θ) = ([θ[1]^2 + θ[2]^2 - 2, θ[1] - θ[2]], [2θ[1] 2θ[2]; 1 -1])
    res = setprecision(BigFloat, 300) do
        CR.gauss_newton(F, BigFloat[1.3, 0.8]; step_tol = big(2.0)^-290, res_floor = big(2.0)^-290,
                        rank_rtol = big(2.0)^-150)
    end
    @test res.converged
    @test abs(res.θ[1] - 1) < big(2.0)^-280
    @test res.iterations <= 10
    @test res.cond ≈ 2 atol = 1e-6            # J = [2 2; 1 -1] at (1, 1): σ = √8, √2
    # rank-deficient (underdetermined) system: pseudoinverse step onto the solution set
    G(θ) = ([θ[1] + θ[2] - 1], [1.0 1.0])
    r2 = CR.gauss_newton(G, [0.0, 0.0]; step_tol = 1e-15, res_floor = 1e-15, rank_rtol = 1e-12)
    @test r2.converged
    @test sum(r2.θ) ≈ 1
end

@testset "rank-revealing solve in extended precision" begin
    J = setprecision(() -> BigFloat[1 2; 2 4; 1 1], 256)
    r = setprecision(() -> BigFloat[3, 6, 2], 256)
    Δ, κ, rank = setprecision(() -> CR.lsq_step(J, r; rank_rtol = big(2.0)^-100), 256)
    @test rank == 2
    @test Float64(norm(J * Δ - r)) < 1e-60
    Js = setprecision(() -> BigFloat[1 2; 2 4], 256)
    _, κs, ranks = setprecision(() -> CR.lsq_step(Js, BigFloat[1, 2]; rank_rtol = big(2.0)^-100), 256)
    @test ranks == 1
    @test κs > 1e20
end

@testset "guard digits from the measured condition number" begin
    @test CR.guard_bits_from_cond(1.0) == 32
    @test CR.guard_bits_from_cond(1e8) == 32 + 32
    @test CR.guard_bits_from_cond(2.0^20) == 32 + 24
    r = rule(XiaoGimbutas(), Simplex{2}(); degree = 15, digits = 60)
    c = certificate(r)
    @test c.cond > 1
    @test c.guard_digits >= floor(Int, CR.guard_bits_from_cond(c.cond) * log10(2)) - 3
    @test c.digits == 60
    @test c.residual < big(10.0)^-59
end

@testset "cancellation token" begin
    tok = CancellationToken()
    @test !CR.iscancelled(tok)
    cancel!(tok)
    @test CR.iscancelled(tok)
    @test_throws CancelledError rule(XiaoGimbutas(), Simplex{2}(); degree = 12, digits = 50, cancel = tok)
    @test_throws CancelledError CR.multistart(CR.xg_entry_for(6).structure, 6; nstarts = 64, cancel = tok)
    # an uncancelled token changes nothing
    @test rule(Simplex{2}(); degree = 12, cancel = CancellationToken()) == rule(Simplex{2}(); degree = 12)
end

@testset "divergence is an error, never a silently unrefined seed" begin
    e = CR.xg_entry_for(10)
    bad = fill(0.02, length(e.seed))      # nowhere near a rule
    @test_throws RefinementError rule(XiaoGimbutas(), Simplex{2}(); degree = 10, seed = ExplicitSeed(bad))
    @test_throws ArgumentError rule(XiaoGimbutas(), Simplex{2}(); degree = 10, seed = ExplicitSeed([0.1]))
    # the shipped seed, passed explicitly, gives the same rule as the table
    r = rule(XiaoGimbutas(), Simplex{2}(); degree = 10, seed = ExplicitSeed(e.seed; source = "test"))
    @test r == rule(XiaoGimbutas(), Simplex{2}(); degree = 10)
    @test occursin("test", provenance(r).seed_source)
end

@testset "seeds: orbit structure alone recovers a rule (licence-free route)" begin
    e = CR.xg_entry_for(6)
    found = CR.multistart(e.structure, 6; nstarts = 256, first_only = false)
    @test !isempty(found)
    @test all(f -> f[2] > 0 && f[3] > 0, found)
    r = rule(XiaoGimbutas(), Simplex{2}(); degree = 6, seed = CR.MultistartSeed(nstarts = 256))
    @test passed(check(r))
    @test occursin("multistart", provenance(r).seed_source)
    @test length(CR.candidate_structures(25, 10)) >= 1
    @test all(s -> CR.npoints(s) == 25, CR.candidate_structures(25, 10))
end

@testset "node elimination and growing from the degree below" begin
    # triangle: grow the degree-7 rule to degree 8, then eliminate
    e7 = CR.xg_entry_for(7)
    B8 = CR.invariant_basis(3, 8)
    g = CR.grow(e7.structure, e7.seed, 8; basis = B8)
    @test g !== nothing
    s, θ = g
    @test CR.npoints(s) > 16
    se, θe = CR.eliminate(s, θ, 8; basis = B8)
    @test CR.npoints(se) < CR.npoints(s)
    @test CR.nunknowns(se) >= CR.n_invariants(3, 8) - 2
    wmin, λmin = CR.rule_margins(se, θe)
    @test wmin > 0 && λmin > 0
    sys = CR.SymmetricMomentSystem(se, 8, Float64, B8)
    @test maximum(abs, sys(θe)[1]) < 1e-12
    # every move keeps the parameter vector consistent with its structure
    for mv in CR.elimination_moves(s, θ)
        @test CR.nunknowns(mv.structure) == length(mv.θ0)
        @test CR.npoints(mv.structure) == CR.npoints(s) - mv.removed
    end
    # the multi-chain driver is deterministic and returns a verifiable rule
    r1 = CR.grow_and_eliminate(e7.structure, e7.seed, 8; chains = 3, basis = B8)
    r2 = CR.grow_and_eliminate(e7.structure, e7.seed, 8; chains = 3, basis = B8)
    @test r1[2] == r2[2]
    # as a seed source
    r = rule(XiaoGimbutas(), Simplex{2}(); degree = 8, seed = LowerDegreeSeed(chains = 3))
    @test passed(check(r))
    @test occursin("grown from the degree below", provenance(r).seed_source)
    @test_throws ArgumentError rule(XiaoGimbutas(), Simplex{2}(); degree = 1, seed = LowerDegreeSeed())
end
