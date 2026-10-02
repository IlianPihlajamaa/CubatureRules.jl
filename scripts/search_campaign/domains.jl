# Domain adapters for the pooled search (see README.md). The engine never touches orbit
# structures directly: everything domain-specific goes through these functions, each a thin
# wrapper over machinery the package already has.
#
#     table_file(d)            shipped seed table, in src/data
#     degree_step(d)           1 on simplices; 2 on boxes, whose rules have odd degree only
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
domain_name(d::SimplexDomain) = d.N == 3 ? "triangle" : "tetrahedron"
table_file(d::SimplexDomain) = d.N == 3 ? "triangle_s3_seeds.toml" : "tetrahedron_s4_seeds.toml"
degree_step(::SimplexDomain) = 1
n_equations(d::SimplexDomain, n) = CR.n_invariants(d.N, n)
structure(d::SimplexDomain, mults) = CR.SymmetricStructure([Int.(m) for m in mults], d.N)
orbit_mults(::SimplexDomain, s) = [o.mult for o in s.orbits]
# general orbit first, centroid last, as the package's own grow orders them
const _SIMPLEX_TYPES = Dict(N => sort(CR.orbit_pattern_types(N); by = p -> (-length(p.mult), p.mult)) for N in (3, 4))
family(d::SimplexDomain, s) = Tuple(count(o -> o.mult == t.mult, s.orbits) for t in _SIMPLEX_TYPES[d.N])
basis(d::SimplexDomain, n) = CR.invariant_basis(d.N, n)

function fit(d::SimplexDomain, s, θ0, n)
    θ = CR.fit_symmetric(s, θ0, n, basis(d, n))
    return θ === nothing ? nothing : (s, θ)
end
system(d::SimplexDomain, s, n) = CR.SymmetricMomentSystem(s, n, Float64, basis(d, n); representatives = true)
weight_indices(::SimplexDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(::SimplexDomain, s, θ) = CR.rule_margins(s, θ)[2] > -0.05
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
margin(::SimplexDomain, s, θ) = CR.rule_margins(s, θ)
function chain(d::SimplexDomain, s, θ, n, seed)
    r = CR.grow_and_eliminate(s, θ, n; basis = basis(d, n), chains = 1, rng_seed = seed, min_excess = 0)
    return r === nothing ? nothing : (r[1], r[2])
end

function refine(d::SimplexDomain, s, n, θ64)
    θ, res, _ = CR.refine_symmetric(s, n, θ64, BITS; basis = basis(d, n))
    return Float64.(θ), res.converged, res.cond
end
function entry(d::SimplexDomain, n, s, θr, cond, note)
    wmin, λmin = CR.rule_margins(s, θr)
    return Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                            "seed" => θr, "status" => "ok", "note" => note, "method" => "pooled search",
                            "min_weight" => wmin, "min_barycentric" => λmin, "cond" => cond)
end

# --- boxes: square (D₄) and cube (O_h) on [-1, 1]^D ---------------------------------------------

struct BoxDomain <: CampaignDomain
    D::Int
end
domain_name(d::BoxDomain) = d.D == 2 ? "square" : "cube"
table_file(d::BoxDomain) = d.D == 2 ? "square_d4_seeds.toml" : "cube_oh_seeds.toml"
degree_step(::BoxDomain) = 2
n_equations(d::BoxDomain, n) = length(CR.box_exponents(d.D, n))
structure(d::BoxDomain, mults) = CR.BoxStructure(d.D, [CR.BoxOrbit(Int.(m), d.D) for m in mults])
orbit_mults(::BoxDomain, s) = [o.mult for o in s.orbits]
const _BOX_TYPES = Dict(D => CR.box_orbit_types(D) for D in (2, 3))
family(d::BoxDomain, s) = Tuple(count(o -> o.mult == t.mult, s.orbits) for t in _BOX_TYPES[d.D])

fit(::BoxDomain, s, θ0, n) = CR.box_fit(s, θ0, n)
system(::BoxDomain, s, n) = CR.BoxMomentSystem(s, n, Float64)
weight_indices(::BoxDomain, s) = [off + 1 for off in CR.param_offsets(s)]
inside(::BoxDomain, s, θ) = all(o_off -> all(abs(θ[o_off[2] + 1 + i]) < 1.05 for i in 1:CR.nparams(o_off[1])),
                                zip(s.orbits, CR.param_offsets(s)))
moves(::BoxDomain, s, θ) = [(mv[1], mv[2], mv[3], mv[4], "box move") for mv in CR.box_moves(s, θ)]

general_orbit(d::BoxDomain) = last(_BOX_TYPES[d.D])
small_orbits(d::BoxDomain) = filter(t -> !isempty(t.mult), _BOX_TYPES[d.D][1:(end - 1)])   # never a second centre
volume(d::BoxDomain) = 2.0^d.D
random_orbit(d::BoxDomain, o, rng) = (o, vcat(0.02 * volume(d) / CR.orbit_size(o), 0.03 .+ 0.94 .* rand(rng, CR.nparams(o))))
function with_orbits(d::BoxDomain, s, θ, new)
    blocks = [(o, θ[(off + 1):(off + CR.nunknowns(o))]) for (o, off) in zip(s.orbits, CR.param_offsets(s))]
    return CR.from_box_blocks(d.D, vcat(blocks, new))
end
margin(::BoxDomain, s, θ) = (m = CR.box_margins(s, θ); (m[1], min(m[2], m[3])))
function chain(::BoxDomain, s, θ, n, seed)
    r = CR.box_grow_and_eliminate(s, θ, n; chains = 1, rng_seed = seed)
    return r === nothing ? nothing : (r[1], r[2])
end

function refine(::BoxDomain, s, n, θ64)
    θ, res, _ = CR.refine_box(s, n, θ64, BITS)
    return Float64.(θ), res.converged, res.cond
end
entry(d::BoxDomain, n, s, θr, cond, note) =
    Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => orbit_mults(d, s),
                     "seed" => θr, "status" => "ok", "note" => note)

const DOMAINS = Dict("triangle" => SimplexDomain(3), "tetrahedron" => SimplexDomain(4),
                     "square" => BoxDomain(2), "cube" => BoxDomain(3))
