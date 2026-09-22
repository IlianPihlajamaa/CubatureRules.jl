using CubatureRules, Test
const CR = CubatureRules

@testset "octahedral elimination moves" begin
    s = CR.OctahedralStructure([:a1, :b, :d])
    θ = [0.05, 0.03, 0.35, 0.04, 0.3, 0.5]
    moves = CR.octahedral_elimination_moves(s, θ)
    @test !isempty(moves)
    @test issorted([-mv.removed for mv in moves])           # largest saving first
    # every move really is smaller, and the parameters match the new structure
    for mv in moves
        @test CR.npoints(mv.structure) == CR.npoints(s) - mv.removed
        @test length(mv.θ0) == CR.nunknowns(mv.structure)
        @test all(>(0), mv.θ0[CR.param_offsets(mv.structure) .+ 1])   # weights stay positive
    end
    @test any(mv -> occursin("drop", mv.what), moves)
    @test any(mv -> occursin("d to b", mv.what), moves)
    # a parameterless orbit is unique, so nothing collapses onto one already present
    @test !any(mv -> occursin("to a1", mv.what), moves)
    s2 = CR.OctahedralStructure([:b])
    @test any(mv -> occursin("to a1", mv.what), CR.octahedral_elimination_moves(s2, [0.05, 0.3]))
    # dropping an orbit preserves the total mass, which is what keeps the fit close
    drop = first(filter(mv -> occursin("drop", mv.what), moves))
    mass(st, p) = sum(p[off + 1] * CR.orbit_size(o) for (o, off) in zip(st.orbits, CR.param_offsets(st)))
    @test mass(drop.structure, drop.θ0) ≈ mass(s, θ) rtol = 1e-12
end

@testset "grow recipes" begin
    e = CR.lebedev_entry_for(15)
    rs = CR.grow_recipes(e.structure, 17)
    @test !isempty(rs)
    m = length(CR.invariant_exponents(17))
    # every recipe reaches at least as many unknowns as there are equations
    for c in rs
        @test CR.nunknowns(e.structure) + sum(CR.nunknowns(CR.OctahedralOrbit(k)) for k in c) >= m
    end
    # cheapest in added points first, and a single 24-point orbit leads
    added(c) = sum(CR.orbit_size(CR.OctahedralOrbit(k)) for k in c)
    @test issorted(added.(rs))
    @test added(first(rs)) == 24
    # a parameterless orbit already present is never offered again
    @test !any(c -> :a1 in c, rs)                 # the degree-15 rule already has a1
end

@testset "grow and eliminate reaches the known degree-17 count" begin
    e = CR.lebedev_entry_for(15)
    got = CR.grow_and_eliminate_octahedral(e.structure, e.seed, 17; chains = 2)
    @test got !== nothing
    s, θ, counts = got
    @test CR.npoints(s) == 110                     # Lebedev's count, reached from the rule below
    @test all(==(110), counts)
    # and it is a real rule once refined, which is what the seed is for: the Float64 fit
    # stops at 1e-12 on the invariants, below what verification at Float64 demands
    bits = CR.digits_to_bits(30)
    θb, res, guard = CR.refine_octahedral(s, 17, θ, bits)
    @test res.converged
    # rounded to the output precision, as the family does through `finalize_number`. The
    # weights are copied out of θ without arithmetic, so without rounding they keep their
    # guard bits and verification holds the rule to an accuracy it never aimed at.
    xs0, ws0 = CR.with_bits(() -> CR.expand(s, θb), bits + guard)
    xs = [CR.SVector{3,BigFloat}(BigFloat(x[1]; precision = bits), BigFloat(x[2]; precision = bits),
                                 BigFloat(x[3]; precision = bits)) for x in xs0]
    ws = [BigFloat(w; precision = bits) for w in ws0]
    r = CR.QuadratureRule(xs, ws, Sphere{3}(), PolynomialDegree(17),
                          CR.Provenance(family = "test", derivation = CR.Derived(), symmetry = :Oh))
    v = check(r)
    @test v.exact && v.sharp === true && v.positive && v.symmetric === true
    @test CR.octahedral_margins(s, θ)[1] > 0
    @test CR.node_separation(s, θ) > 1e-6
end

@testset "fit_octahedral" begin
    e = CR.lebedev_entry_for(17)
    # refitting a solution at its own degree returns it
    θ = CR.fit_octahedral(e.structure, e.seed, 17)
    @test θ !== nothing
    @test maximum(abs, θ - CR.canonical_octahedral(e.structure, e.seed)) < 1e-9
    # a structure that cannot meet the degree is rejected rather than half-fitted
    @test CR.fit_octahedral(CR.OctahedralStructure([:a1]), [0.5], 9) === nothing
end
