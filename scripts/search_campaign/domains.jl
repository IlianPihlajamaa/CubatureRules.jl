# Domain adapters for the pooled search (see README.md). The engine never touches orbit
# structures directly: everything domain-specific goes through these functions, each a thin
# wrapper over machinery the package already has.
#
#     table_file(d)            shipped seed table, in src/data
#     degree_step(d)           1 on simplices; 2 on boxes and the sphere, whose rules have odd
#                              degree only
#     n_equations(d, n)        invariant moment equations at degree n
#     structure(d, mults)      orbit structure from the table's `structure` field, and back
#     orbit_mults(d, s)
#     family(d, s)             orbit-type counts: the key the finish step balances over
#     fit(d, s, θ0, n)         Float64 fit onto the degree-n moment variety, accepted only as a
#                              positive interior rule: canonical (s, θ), or nothing
#     system(d, s, n)          Float64 moment system (r, J) = sys(θ), for the positive-weight solve
#     weight_indices(d, s)     positions of the weights in θ
#     inside(d, s, θ)          loose admissibility test used during a solve
#     moves(d, s, θ)           elimination moves: (s, θ0, removed points, score, label)
#     general_orbit(d), small_orbits(d), random_orbit(d, o, rng), volume(d)   for growing
#     with_orbits(d, s, θ, new)  (s, θ) with orbit blocks `new` added, canonical
#     margin(d, s, θ)          (smallest weight, smallest interior margin)
#     chain(d, s, θ, n, seed)  one of the package's own grow → eliminate chains: (s, θ) or nothing
#     refine(d, s, n, θ64)     100-digit refinement: (θ refined, converged, cond)
#     entry(d, n, s, θr, cond, note)   a table row in the shipped format

const BITS = CR.digits_to_bits(100)

abstract type CampaignDomain end

# --- simplices: triangle (S₃) and tetrahedron (S₄) -----------------------------------------

struct SimplexDomain <: CampaignDomain
    N::Int
end
domain_name(d::SimplexDomain) = ("triangle", "tetrahedron", "simplex4")[d.N - 2]
table_file(d::SimplexDomain) = ("triangle_s3_seeds.toml", "tetrahedron_s4_seeds.toml", "simplex4_s5_seeds.toml")[d.N - 2]
degree_step(::SimplexDomain) = 1
n_equations(d::SimplexDomain, n) = CR.n_invariants(d.N, n)
structure(d::SimplexDomain, mults) = CR.SymmetricStructure([Int.(m) for m in mults], d.N)
orbit_mults(::SimplexDomain, s) = [o.mult for o in s.orbits]
# general orbit first, centroid last, as the package's own grow orders them
const _SIMPLEX_TYPES = Dict(N => sort(CR.orbit_pattern_types(N); by = p -> (-length(p.mult), p.mult)) for N in (3, 4, 5))
family(d::SimplexDomain, s) = Tuple(count(o -> o.mult == t.mult, s.orbits) for t in _SIMPLEX_TYPES[d.N])
basis(d::SimplexDomain, n) = CR.invariant_basis(d.N, n)

function fit(d::SimplexDomain, s, θ0, n)
    θ = CR.fit_symmetric(s, θ0, n, basis(d, n))
    return θ === nothing ? nothing : (s, θ)
end
system(d::SimplexDomain, s, n) = CR.SymmetricMomentSystem(s, n, Float64, basis(d, n); representatives = true)
weight_indices(::SimplexDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(::SimplexDomain, s, θ) = CR.margins(s, θ).interior > -0.05
moves(::SimplexDomain, s, θ) = [(mv.structure, mv.θ0, mv.removed, mv.score, mv.what) for mv in CR.elimination_moves(s, θ)]

general_orbit(d::SimplexDomain) = _SIMPLEX_TYPES[d.N][1]
small_orbits(d::SimplexDomain) = _SIMPLEX_TYPES[d.N][2:(end - 1)]          # never a second centroid
volume(d::SimplexDomain) = 1 / factorial(d.N - 1)
function random_orbit(d::SimplexDomain, p, rng)
    u = -log.(rand(rng, length(p.mult)))
    u ./= sum(u)
    return (copy(p.mult), vcat(0.02 * volume(d) / CR.orbit_size(p), [u[i] / p.mult[i] for i in 1:(length(p.mult) - 1)]))
end
with_orbits(d::SimplexDomain, s, θ, new) = CR.from_blocks(d.N, vcat(CR.orbit_blocks(s, θ), new))
margin(::SimplexDomain, s, θ) = CR.margins(s, θ)
function chain(d::SimplexDomain, s, θ, n, seed)
    r = CR.grow_and_eliminate(s, θ, n; basis = basis(d, n), chains = 1, rng_seed = seed, min_excess = 0)
    return r === nothing ? nothing : (r[1], r[2])
end

function refine(d::SimplexDomain, s, n, θ64)
    θ, res, _ = CR.refine_seed(s, n, θ64, BITS; basis = basis(d, n))
    return Float64.(θ), res.converged, res.cond
end
function entry(d::SimplexDomain, n, s, θr, cond, note)
    wmin, λmin = CR.margins(s, θr)
    return Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                            "seed" => θr, "status" => "ok", "note" => note, "method" => "pooled search",
                            "min_weight" => wmin, "min_barycentric" => λmin, "cond" => cond)
end

# --- boxes and balls: hyperoctahedral orbits ([-1, 1]^D for D = 2, 3, 4; the disk and the ball) ----
#
# The same orbits and moves on the box and on the ball; the package switches the moment system,
# margins and random values by a geometry (refine/box_search.jl).

struct BoxDomain <: CampaignDomain
    D::Int
    round::Bool
end
BoxDomain(D) = BoxDomain(D, false)
geom(d::BoxDomain) = d.round ? CR.RoundGeometry() : CR.BoxGeometry()
domain_name(d::BoxDomain) = d.round ? ("disk", "ball")[d.D - 1] : ("square", "cube", "hypercube4")[d.D - 1]
table_file(d::BoxDomain) = d.round ? ("disk_d4_seeds.toml", "ball_oh_seeds.toml")[d.D - 1] :
                           ("square_d4_seeds.toml", "cube_oh_seeds.toml", "hypercube_b4_seeds.toml")[d.D - 1]
degree_step(::BoxDomain) = 2
n_equations(d::BoxDomain, n) = CR.orbit_equations(geom(d), d.D, n)
structure(d::BoxDomain, mults) = CR.BoxStructure(d.D, [CR.BoxOrbit(Int.(m), d.D) for m in mults])
orbit_mults(::BoxDomain, s) = [o.mult for o in s.orbits]
const _BOX_TYPES = Dict(D => CR.box_orbit_types(D) for D in (2, 3, 4))
family(d::BoxDomain, s) = Tuple(count(o -> o.mult == t.mult, s.orbits) for t in _BOX_TYPES[d.D])

fit(d::BoxDomain, s, θ0, n) = CR.box_fit(s, θ0, n; geom = geom(d))
system(d::BoxDomain, s, n) = CR.orbit_system(geom(d), s, n, Float64)
weight_indices(::BoxDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(d::BoxDomain, s, θ) = all(o_off -> CR.loosely_inside(geom(d), o_off[1], [θ[o_off[2] + 1 + i] for i in 1:CR.nparams(o_off[1])]),
                                 zip(s.orbits, CR.param_offsets(s)))
moves(::BoxDomain, s, θ) = [(mv[1], mv[2], mv[3], mv[4], "box move") for mv in CR.box_moves(s, θ)]

general_orbit(d::BoxDomain) = last(_BOX_TYPES[d.D])
small_orbits(d::BoxDomain) = filter(t -> !isempty(t.mult), _BOX_TYPES[d.D][1:(end - 1)])   # never a second centre
volume(d::BoxDomain) = CR.orbit_volume(geom(d), d.D)
random_orbit(d::BoxDomain, o, rng) = (o, vcat(0.02 * volume(d) / CR.orbit_size(o), CR.random_values(geom(d), o.mult, rng)))
function with_orbits(d::BoxDomain, s, θ, new)
    blocks = [(o, θ[(off + 1):(off + CR.nunknowns(o))]) for (o, off) in zip(s.orbits, CR.param_offsets(s))]
    return CR.from_box_blocks(d.D, vcat(blocks, new))
end
margin(::BoxDomain, s, θ) = (m = CR.margins(s, θ); (m.weight, min(m.interior, m.gap)))
function chain(::BoxDomain, s, θ, n, seed)
    r = CR.box_grow_and_eliminate(s, θ, n; chains = 1, rng_seed = seed)
    return r === nothing ? nothing : (r[1], r[2])
end

function refine(::BoxDomain, s, n, θ64)
    θ, res, _ = CR.refine_seed(s, n, θ64, BITS)
    return Float64.(θ), res.converged, res.cond
end
entry(d::BoxDomain, n, s, θr, cond, note) =
    Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                     "seed" => θr, "status" => "ok", "note" => note)

# --- the pyramid (C₄ᵥ) and the wedge (D₃ₕ), through the package's generic orbit search ---------
#
# refine/orbit_search.jl defines one interface over PyramidStructure and WedgeStructure
# (orbit_types, orbit_fit, structure_moves, orbit_grow_and_eliminate, …); this is a thin layer
# over it.

struct OrbitDomain{T} <: CampaignDomain end
domain_name(::OrbitDomain{CR.PyramidStructure}) = "pyramid"
domain_name(::OrbitDomain{CR.WedgeStructure}) = "wedge"
table_file(::OrbitDomain{CR.PyramidStructure}) = "pyramid_c4v_seeds.toml"
table_file(::OrbitDomain{CR.WedgeStructure}) = "wedge_d3h_seeds.toml"
degree_step(::OrbitDomain) = 1
n_equations(::OrbitDomain{T}, n) where {T} = CR.n_moment_equations(T, n)
structure(::OrbitDomain{T}, lists) where {T} = T([Int.(l) for l in lists])
orbit_mults(::OrbitDomain{CR.PyramidStructure}, s) = [o.base.mult for o in s.orbits]
orbit_mults(::OrbitDomain{CR.WedgeStructure}, s) = CR.wedge_structure_lists(s)
family(::OrbitDomain{T}, s) where {T} = Tuple(count(==(t), s.orbits) for t in CR.orbit_types(T))

fit(::OrbitDomain, s, θ0, n) = CR.orbit_fit(s, θ0, n)
system(::OrbitDomain, s, n) = CR.orbit_system(s, n, Float64)
weight_indices(::OrbitDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(::OrbitDomain, s, θ) = all(o_off -> CR.loosely_inside(o_off[1], θ[(o_off[2] + 2):(o_off[2] + CR.nunknowns(o_off[1]))]),
                                  zip(s.orbits, CR.param_offsets(s)))
moves(::OrbitDomain, s, θ) = [(mv[1], mv[2], mv[3], mv[4], "orbit move") for mv in CR.structure_moves(s, θ)]

general_orbit(::OrbitDomain{T}) where {T} = CR.general_orbit(T)
# orbits with a free parameter can repeat; one without (the wedge's centroid) is unique
small_orbits(d::OrbitDomain{T}) where {T} = filter(o -> o != general_orbit(d) && CR.nparams(o) > 0, CR.orbit_types(T))
volume(::OrbitDomain{T}) where {T} = CR.structure_volume(T)
random_orbit(d::OrbitDomain, o, rng) = (o, vcat(0.02 * volume(d) / CR.orbit_size(o), CR.random_params(o, rng)))
function with_orbits(::OrbitDomain{T}, s, θ, new) where {T}
    blocks = vcat([(o, θ[(off + 1):(off + CR.nunknowns(o))]) for (o, off) in zip(s.orbits, CR.param_offsets(s))], new)
    return CR.canonicalize(CR.make_structure(T, [b[1] for b in blocks]), reduce(vcat, (b[2] for b in blocks); init = Float64[]))
end
margin(::OrbitDomain, s, θ) = (m = CR.structure_margins(s, θ); (m[1], min(m[2], m[3])))
function chain(::OrbitDomain, s, θ, n, seed)
    r = CR.orbit_grow_and_eliminate(s, θ, n; chains = 1, rng_seed = seed)
    return r === nothing ? nothing : (r[1], r[2])
end

function refine(::OrbitDomain{T}, s, n, θ64) where {T}
    θ, res, _ = (T === CR.PyramidStructure ? CR.refine_pyramid : CR.refine_wedge)(s, n, θ64, BITS)
    return Float64.(θ), res.converged, res.cond
end
entry(d::OrbitDomain, n, s, θr, cond, note) =
    Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                     "seed" => θr, "status" => "ok", "note" => note)

# --- the sphere: O_h-symmetric (Lebedev-type) rules ------------------------------------------

struct SphereDomain <: CampaignDomain end

const _OCT_KINDS = (:a1, :a2, :a3, :b, :c, :d)

table_file(::SphereDomain) = "lebedev_seeds.toml"
degree_step(::SphereDomain) = 2
n_equations(::SphereDomain, n) = length(CR.invariant_exponents(n))
structure(::SphereDomain, kinds) = CR.OctahedralStructure(Symbol.(kinds))
orbit_mults(::SphereDomain, s) = [String(o.kind) for o in s.orbits]
family(::SphereDomain, s) = Tuple(count(o -> o.kind === k, s.orbits) for k in _OCT_KINDS)

# The fit works in the orthonormal invariant harmonics, not in p₄ᵃp₆ᵇ as the package's
# fit_octahedral does: in Float64 the monomial system is hopeless above degree 17 (condition
# number 6e5 at the degree-17 rule, 1e19 one orbit later, 1e48 at degree 133), the harmonic
# one is not (23 and 2.6e3).
# A b orbit needs l < 1/√2 and a d orbit r² + s² < 1; a solve or a perturbation can step past
# that, where the system throws a DomainError. Answered with a huge residual instead, a solve
# steps back rather than failing the whole task.
struct SafeSphereSystem{T}
    sys::T
end
function (w::SafeSphereSystem)(θ; jacobian::Bool = true)
    try
        return w.sys(θ; jacobian)
    catch e
        e isa DomainError || rethrow()
        m, u = CR.n_equations(w.sys), CR.n_unknowns(w.sys)
        return fill(1e10, m), jacobian ? zeros(m, u) : zeros(0, 0)
    end
end
sphere_inside(s, θ) = try CR.octahedral_margins(s, θ)[2] > -0.02 catch e; e isa DomainError || rethrow(); false end

function fit(::SphereDomain, s, θ0, n)
    sphere_inside(s, θ0) || return nothing             # a start outside an orbit's domain is no start
    sys = SafeSphereSystem(CR.OctahedralHarmonicSystem(s, n, Float64))
    θ, nr = CR.levenberg_marquardt(sys, θ0; maxiter = 400, tol = 1e-13, accept = x -> sphere_inside(s, x))
    nr <= 1e-10 || return nothing
    res = CR.gauss_newton(sys, θ; step_tol = 1e-15, res_floor = 1e-13, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    wmin, dmin = CR.octahedral_margins(s, res.θ)
    (wmin > 0 && dmin > 1e-6 && CR.node_separation(s, res.θ) > 1e-6) || return nothing
    return s, CR.canonical_octahedral(s, res.θ)
end
system(::SphereDomain, s, n) = SafeSphereSystem(CR.OctahedralHarmonicSystem(s, n, Float64))
weight_indices(::SphereDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(::SphereDomain, s, θ) = sphere_inside(s, θ)
moves(::SphereDomain, s, θ) = [(mv.structure, mv.θ0, mv.removed, 0.0, mv.what) for mv in CR.octahedral_elimination_moves(s, θ)]

general_orbit(::SphereDomain) = CR.OctahedralOrbit(:d)
small_orbits(::SphereDomain) = [CR.OctahedralOrbit(:b), CR.OctahedralOrbit(:c)]   # a₁, a₂, a₃ occur at most once
volume(::SphereDomain) = 4π
random_orbit(::SphereDomain, o, rng) = (o.kind, CR.random_octahedral_parameters(rng, CR.OctahedralStructure([o.kind]), 1.0))
# New orbits start at half the rule's mean weight, as in the package's grow_octahedral: sphere
# degrees step by two, and an orbit carrying a token weight cannot absorb the gap.
function with_orbits(::SphereDomain, s, θ, new)
    w = 0.5 * 4π / CR.npoints(s)
    return CR.from_octahedral_blocks(vcat(CR.octahedral_blocks(s, θ), [(k, vcat(w, p[2:end])) for (k, p) in new]))
end
margin(::SphereDomain, s, θ) = CR.octahedral_margins(s, θ)
# A grow → eliminate chain with the fit above (the package's grow_and_eliminate_octahedral uses
# the monomial fit). Growth always adds a d orbit: from degree 18 on every exact rule needs one,
# since the square of the product of the nine mirror planes' linear forms is an invariant of
# degree 18 that vanishes, with its gradient, at every point of every other orbit type, and has
# a positive integral.
function chain(d::SphereDomain, s, θ, n, seed)
    rng = MersenneTwister(seed)
    m = n_equations(d, n)
    g = nothing
    for _ in 1:16
        new = [random_orbit(d, general_orbit(d), rng)]
        u = CR.nunknowns(s) + CR.nunknowns(general_orbit(d))
        while u < m
            o = rand(rng, small_orbits(d))
            push!(new, random_orbit(d, o, rng)); u += CR.nunknowns(o)
        end
        g = fit(d, with_orbits(d, s, θ, new)..., n)
        g === nothing || break
    end
    g === nothing && return nothing
    s, θ = g
    while true                           # eliminate: most points first, random among equals
        mvs = sort!(filter(mv -> CR.nunknowns(mv[1]) >= m, moves(d, s, θ)); by = mv -> (-mv[3], rand(rng)))
        r = nothing
        for mv in mvs
            r = fit(d, mv[1], mv[2], n)
            r === nothing || break
        end
        r === nothing && return s, θ
        s, θ = r
    end
end

function refine(::SphereDomain, s, n, θ64)
    θ, res, _ = CR.refine_octahedral(s, n, θ64, BITS)
    return Float64.(θ), res.converged, res.cond
end
function entry(d::SphereDomain, n, s, θr, cond, note)
    wmin, dmin = CR.octahedral_margins(s, θr)
    return Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                            "seed" => θr, "status" => "ok", "note" => note, "method" => "pooled search",
                            "min_weight" => wmin, "min_margin" => dmin, "cond" => cond)
end

const DOMAINS = Dict("triangle" => SimplexDomain(3), "tetrahedron" => SimplexDomain(4), "simplex4" => SimplexDomain(5),
                     "square" => BoxDomain(2), "cube" => BoxDomain(3), "hypercube4" => BoxDomain(4),
                     "disk" => BoxDomain(2, true), "ball" => BoxDomain(3, true),
                     "pyramid" => OrbitDomain{CR.PyramidStructure}(), "wedge" => OrbitDomain{CR.WedgeStructure}(),
                     "sphere" => SphereDomain())
