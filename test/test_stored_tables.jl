using CubatureRules, Test
const CR = CubatureRules

# Float64 requests for the seeded families are served from the stored tables without
# refinement, so the tables themselves must be right. scripts/certify_tables.jl records each
# entry's residual; here every recorded residual is recomputed from scratch — in the full
# sweep (CUBATURERULES_FULL_SWEEP=1) for every entry, and per commit for the simplex and box
# entries of at most 650 points, as in the other table sweeps. Above that the tetrahedron
# entries need invariant bases that take a minute each to build (70 s at degree 30).
stored_per_commit(e) = get(ENV, "CUBATURERULES_FULL_SWEEP", "0") == "1" || e.npoints <= 650

@testset "every stored seed has a recorded residual, and it is reproduced" begin
    tables = [(3, filter(stored_per_commit, CR.xg_entries())), (4, filter(stored_per_commit, CR.tet_entries())),
              (0, CR.lebedev_entries())]
    for (N, entries) in tables, e in entries
        @test isfinite(e.residual) && e.residual < 1e-14 && e.residual_bits >= 256
        r = CR.with_bits(128) do
            θ = BigFloat.(e.seed)
            res = N == 0 ? CR.OctahedralMomentSystem(e.structure, e.degree, BigFloat)(θ)[1] :
                  CR.SymmetricMomentSystem(e.structure, e.degree, BigFloat, CR.invariant_basis(N, e.degree))(θ; jacobian = false)[1]
            Float64(maximum(abs, res))
        end
        @test isapprox(r, e.residual; rtol = 1e-6, atol = 1e-30)
    end
    for D in (2, 3), e in filter(stored_per_commit, CR.box_entries(D))
        @test isfinite(e.residual) && e.residual < 1e-14 && e.residual_bits >= 256
        r = CR.with_bits(128) do
            Float64(maximum(abs, CR.BoxMomentSystem(e.structure, e.degree, BigFloat)(BigFloat.(e.seed); jacobian = false)[1]))
        end
        @test isapprox(r, e.residual; rtol = 1e-6, atol = 1e-30)
    end
end

@testset "Float64 is shipped, lower precisions are rounded, higher ones refined" begin
    r = rule(Simplex{2}(); degree = 30)
    c = certificate(r)
    @test c.iterations == 0 && c.guard_digits == 0 && c.residual_bits == 256 && c.residual < 1e-14
    @test occursin("shipped as stored", last(provenance(r).path))
    @test CR.passed(verify(r))
    @test (@elapsed rule(Simplex{2}(); degree = 30)) < 0.5

    # the shipped rule is the refined rule to rounding
    rr = rule(Simplex{2}(); degree = 30, digits = 30)
    @test certificate(rr).iterations >= 1
    @test maximum(maximum(abs.(Float64.(a) .- b)) for (a, b) in zip(nodes(rr), nodes(r))) <= 2eps()
    @test maximum(abs.(Float64.(weights(rr)) .- weights(r)) ./ weights(r)) <= 4eps()

    # below Float64: rounded, with the residual of the rounded parameters evaluated in Float64
    for T in (Float32, Float16)
        rl = rule(Simplex{2}(); degree = 10, T)
        @test eltype(rl) == T
        @test certificate(rl).residual_bits == 53 && 0 < certificate(rl).residual < 10eps(T)
        @test occursin("rounded to $T", last(provenance(rl).path))
    end
    @test CR.passed(verify(rule(Simplex{2}(); degree = 10, T = Float32)))

    # tetrahedra and spheres take the same path
    for dom in (Simplex{3}(), Sphere{3}())
        rs = rule(dom; degree = 9)
        @test certificate(rs).iterations == 0 && CR.passed(verify(rs))
    end

    # any other seed is refined even in Float64
    e = CR.xg_entry_for(10)
    rx = rule(XiaoGimbutas(), Simplex{2}(); degree = 10, seed = CR.ExplicitSeed(e.seed))
    @test any(s -> startswith(s, "refine:"), provenance(rx).path)
end
