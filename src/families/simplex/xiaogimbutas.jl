# Minimal-point fully symmetric triangle rules, refined to arbitrary precision
# (PLAN §6 Tier 3 — the v0.1 cornerstone).
#
# Seeds. No numbers are copied from any published table (PLAN §0.2). The point counts and
# the fact that the rules are fully S_3-symmetric come from Xiao & Gimbutas (2010); the orbit
# structures and Float64 seeds in `src/data/triangle_s3_seeds.toml` were found in-house by
# multistart Gauss–Newton from the orbit structure alone (`scripts/generate_triangle_seeds.jl`)
# and are MIT-licensed. Where several rules share the point count, the one with the largest
# minimum barycentric coordinate is shipped, so a shipped rule need not coincide node for
# node with the published table.

"""
    XiaoGimbutas()

Fully symmetric, positive-weight, interior-node triangle rules with the minimal point counts
of Xiao & Gimbutas (2010), refined by Gauss–Newton on the `S_3`-invariant moment system to
any requested precision. Seeded: available at the degrees in the shipped seed table.

Keyword `seed` of [`rule`](@ref) selects the seed source: `TableSeed()` (default) or
`MultistartSeed()` to recover the rule from its orbit structure alone.
"""
struct XiaoGimbutas <: RuleFamily end

derivation(::Type{XiaoGimbutas}) = Seeded()

const XIAO_GIMBUTAS_2010 = Citation(
    key = "XiaoGimbutas2010", authors = ["Hong Xiao", "Zydrunas Gimbutas"],
    title = "A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions",
    journal = "Computers & Mathematics with Applications", year = 2010, volume = "59", pages = "663--676",
    doi = "10.1016/j.camwa.2009.10.027")

# --- the seed table ---------------------------------------------------------------------

const TRIANGLE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "triangle_s3_seeds.toml")
isfile(TRIANGLE_SEED_FILE) && include_dependency(TRIANGLE_SEED_FILE)

"The fully symmetric seed table of the `(N - 1)`-simplex at `path`."
load_simplex_seeds(path, N) = load_seeds(path, SymmetricStructure, m -> SymmetricStructure(int_lists(m), N))

# Static package data, parsed once at precompile time (not a cache of rules: see §7).
const TRIANGLE_SEEDS = load_simplex_seeds(TRIANGLE_SEED_FILE, 3)

xg_entries() = filter(e -> e.status == "ok", TRIANGLE_SEEDS)
xg_entry_for(degree::Integer) = seed_entry_for(xg_entries(), degree)

function candidates(::Type{XiaoGimbutas}, dom::Simplex{2}, c::PolynomialDegree)
    (isreference(dom) && xg_entry_for(c.d) !== nothing) || return XiaoGimbutas[]
    return [XiaoGimbutas()]
end

npoints(::XiaoGimbutas, dom, degree::Integer) = _xg_entry(degree).npoints
claimed_degree(::XiaoGimbutas, dom, degree) = _xg_entry(degree).degree
function degree_range(::XiaoGimbutas, dom::Simplex{2})
    es = xg_entries()
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::XiaoGimbutas, dom) = 1:0
properties(::XiaoGimbutas, dom, degree) = (positive = true, interior = true, symmetry = SymmetricGroup{3}(), nested = false)
cost_estimate(f::XiaoGimbutas, dom, degree, T) =
    float(npoints(f, dom, degree)) * dubiner_length(claimed_degree(f, dom, degree)) / 10 * _precision_factor(T)

function _xg_entry(degree)
    e = xg_entry_for(degree)
    e === nothing && throw(ArgumentError("no XiaoGimbutas seed covers degree $degree " *
                                         "(shipped range $(degree_range(XiaoGimbutas(), Simplex{2}())))"))
    return e
end

# --- construction --------------------------------------------------------------------------

"What the certificate of a fully symmetric rule on the `(N - 1)`-simplex quotes the residual of."
simplex_equations(N, n) = "$(("S₃", "S₄", "S₅")[N - 2])-invariant moment system in orbit parameters " *
                          "(orthonormal Dubiner basis, degree $n)"
const SIMPLEX_ORIGIN = "generated in-house by multistart from the orbit structure"

function build(f::XiaoGimbutas, dom::Simplex{2}, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed())
    e = _xg_entry(degree)
    return build_seeded(e, ctx; seed, lower = seed_entry_below(xg_entries(), e.degree), family = "XiaoGimbutas",
                        domain = Simplex{2}(), symmetry = SymmetricGroup{3}(), table = "src/data/triangle_s3_seeds.toml",
                        origin = SIMPLEX_ORIGIN, equations = simplex_equations(3, e.degree),
                        citations = [XIAO_GIMBUTAS_2010],
                        license = "MIT (seeds generated in-house; point counts from Xiao & Gimbutas 2010)")
end
