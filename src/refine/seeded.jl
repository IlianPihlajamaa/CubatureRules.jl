# The seed tables behind every fully symmetric family, and seed → refine → certify on them
# (`build_seeded`). A table row holds an orbit structure and Float64 parameters; what differs
# between the domains is the structure, and it answers for itself through `moment_system`,
# `expand` and `margins`.

"""
    SeedEntry{S}

One row of a seed table: the degree it serves, its point count, its orbit structure (an
`S <: OrbitStructure`) and the `Float64` parameters found for it. `cond`, `residual` and
`residual_bits` are recorded by scripts/certify_tables.jl (NaN, NaN and 0 until then);
`method` and `next_error` only where the table records them (the sphere's).
"""
struct SeedEntry{S<:OrbitStructure}
    degree::Int
    npoints::Int
    structure::S
    seed::Vector{Float64}
    status::String
    note::String
    cond::Float64
    residual::Float64            # the moment-system residual of the stored seed, at residual_bits
    residual_bits::Int
    method::String               # how the seed was found
    next_error::Float64          # the miss at degree + 1, measured at high precision
end

"A table row's orbits as lists of multiplicities, each a `Vector{Int}` (TOML stores `[]` untyped)."
int_lists(m) = [Int[x for x in l] for l in m]

"""
    load_seeds(path, S, structure) -> Vector{SeedEntry{S}}

The rows of the TOML seed table at `path`, by degree; `structure` makes an `S` from a row's
`structure` field as stored. Empty if there is no such file.
"""
function load_seeds(path, ::Type{S}, structure) where {S<:OrbitStructure}
    out = SeedEntry{S}[]
    isfile(path) || return out
    for e in get(TOML.parsefile(path), "rule", Any[])
        push!(out, SeedEntry{S}(e["degree"], e["npoints"], structure(e["structure"]), Float64.(e["seed"]),
                                get(e, "status", "ok"), get(e, "note", ""), Float64(get(e, "cond", NaN)),
                                Float64(get(e, "residual", NaN)), Int(get(e, "residual_bits", 0)),
                                get(e, "method", ""), Float64(get(e, "next_error", NaN))))
    end
    return sort!(out; by = e -> e.degree)
end

"The entry serving a request for `degree`: fewest points, then highest degree."
function seed_entry_for(entries, degree::Integer)
    es = filter(e -> e.degree >= degree, entries)
    isempty(es) && return nothing
    return first(sort(es; by = e -> (e.npoints, -e.degree)))
end

"The entry of degree exactly `n - 1`, if shipped."
function seed_entry_below(entries, n::Integer)
    i = findfirst(e -> e.degree == n - 1, entries)
    return i === nothing ? nothing : entries[i]
end

# --- refinement ---------------------------------------------------------------------------

"""
    refine_seed(s, n, θ64, bits; basis, cancel, verbose) -> (θ, result, guard_bits)

Refine a `Float64` seed `θ64` of the fully symmetric rule with orbit structure `s` to `bits`
bits, by [`refine_system`](@ref) on its degree-`n` [`refinement_system`](@ref). `basis` is
what the systems in each arithmetic share ([`moment_basis`](@ref)); pass it to reuse it.
"""
refine_seed(s::OrbitStructure, n::Integer, θ64::Vector{Float64}, bits::Integer; basis = moment_basis(s, n),
            cancel = nothing, verbose::Integer = 0) =
    refine_system(S -> refinement_system(s, n, S, basis), θ64, bits; cancel, verbose,
                  label = "refinement to degree $n ($(nameof(typeof(s))), $(length(s.orbits)) orbits)")

"The largest equation of the degree-`n` [`moment_system`](@ref) of `s` at `θ`, in the arithmetic of `θ`."
moment_residual(s::OrbitStructure, n::Integer, θ::AbstractVector, basis = moment_basis(s, n)) =
    maximum(abs, moment_system(s, n, eltype(θ), basis)(θ; jacobian = false)[1])

"The Cartesian coordinates of a node as `expand` gives it: barycentric on the simplex."
cartesian_node(::OrbitStructure, x) = x
cartesian_node(::SymmetricStructure, λ) = view(λ, 2:length(λ))

"""
    seed_parameters(seed, e, ctx; family, domain, lower, basis, seed_desc) -> (structure, θ64, description)

Where the refinement of row `e` starts: the stored parameters, an explicit vector, a
multistart on the stored structure, or the row one degree down grown and node-eliminated.
"""
seed_parameters(::TableSeed, e::SeedEntry, ctx; seed_desc, kw...) = (e.structure, e.seed, seed_desc)
function seed_parameters(seed::ExplicitSeed, e::SeedEntry, ctx; kw...)
    length(seed.θ) == nunknowns(e.structure) || throw(ArgumentError(
        "an explicit seed for degree $(e.degree) needs $(nunknowns(e.structure)) parameters (structure $(e.structure))"))
    return e.structure, seed.θ, describe(seed)
end
function seed_parameters(seed::MultistartSeed, e::SeedEntry{SymmetricStructure}, ctx; family, basis, kw...)
    found = multistart(e.structure, e.degree; nstarts = seed.nstarts, rng_seed = seed.rng_seed, basis,
                       first_only = false, cancel = ctx.cancel)
    isempty(found) && throw(RefinementError(family,
        "multistart found no valid degree-$(e.degree) rule with structure $(e.structure) in $(seed.nstarts) starts"))
    best = argmax(t -> (t[3], -t[4]), found)
    return e.structure, best[1], describe(seed) * "; $(length(found)) valid rule(s) found, largest minimum barycentric kept"
end
function seed_parameters(seed::MultistartSeed, e::SeedEntry{OctahedralStructure}, ctx; family, kw...)
    found = multistart_octahedral(e.structure, e.degree; nstarts = seed.nstarts, rng_seed = seed.rng_seed,
                                  first_only = false, cancel = ctx.cancel)
    isempty(found) && throw(RefinementError(family,
        "multistart found no valid degree-$(e.degree) rule with structure $(e.structure) in $(seed.nstarts) starts"))
    best = argmax(t -> (t[2], t[3]), found)
    return e.structure, best[1], describe(seed) * "; $(length(found)) valid rule(s) found, largest minimum weight kept"
end
function seed_parameters(seed::LowerDegreeSeed, e::SeedEntry{SymmetricStructure}, ctx; family, lower, basis, kw...)
    n = e.degree
    lower === nothing && throw(ArgumentError("$family has no rule one degree below $n to grow from"))
    got = grow_and_eliminate(lower.structure, lower.seed, n; chains = seed.chains, basis, rng_seed = seed.rng_seed,
                             cancel = ctx.cancel)
    got === nothing && throw(RefinementError(family, "could not grow the degree-$(lower.degree) rule to degree $n"))
    return got[1], got[2], describe(seed) * " from the $(npoints(lower.structure))-point degree-$(lower.degree) rule"
end
seed_parameters(seed::SeedSource, e::SeedEntry, ctx; family, domain, kw...) = throw(ArgumentError(
    "$family on $domain starts from the stored table or an explicit seed, not from $(describe(seed))"))

# --- seed → refine → certify -----------------------------------------------------------------

"""
    build_seeded(e::SeedEntry, ctx; family, domain, symmetry, table, origin, equations,
                 citations = Citation[], comparisons = Citation[], license, seed = TableSeed(),
                 lower = nothing)

Seed → refine → certify for the seed-table row `e`, shared by every family backed by an
orbit-structure table. `table` and `origin` say where the row comes from, `equations` what
the certificate's residual is of. The rule's citations start with this package, whose
searches produced the table; `citations` adds the method it follows, and `comparisons` the
published rules it was compared with. `seed` chooses where the refinement starts (see
[`SeedSource`](@ref)); `lower`, the row one degree down, is what [`LowerDegreeSeed`](@ref)
grows from.

A stored row is a converged `Float64` rule whose residual was checked at high precision when
the table was written. At 53 bits or fewer Newton has nothing to add, so the rule is shipped
as stored, or rounded: its parameters are expanded to nodes and weights at 128 bits and
rounded once. That comes before anything else is computed, since even the simplex's
invariant basis takes seconds at high degree. Otherwise the seed is refined by Gauss–Newton
with guard bits, and the certificate quotes the residual of the delivered parameters at
twice the precision.
"""
function build_seeded(e::SeedEntry, ctx::BuildContext{T}; family::String, domain::Domain, symmetry::SymmetryGroup,
                      table::String, origin::String, equations::String, citations::Vector{Citation} = Citation[],
                      comparisons::Vector{Citation} = Citation[], license::String, seed::SeedSource = TableSeed(),
                      lower::Union{Nothing,SeedEntry} = nothing) where {T}
    isexact(ctx) && throw(ArgumentError("$family nodes are irrational; $(T) is not supported"))
    n = e.degree
    seed_desc = "stored Float64 table $table ($origin)"
    if seed isa TableSeed && ctx.bits <= 53 && isfinite(e.residual)
        ctx.verbose >= 1 && @info "shipping the stored Float64 rule for degree $n (no refinement needed at $(ctx.bits) bits)"
        structure = e.structure
        xs, ws = with_bits(() -> expand(structure, BigFloat.(e.seed)), 128)
        if ctx.bits >= 53
            resid, rbits = e.residual, e.residual_bits
            steps = [@sprintf("shipped as stored, without refinement: residual %.1e at %d bits, checked when the table was written",
                              resid, rbits)]
        else
            resid, rbits = moment_residual(structure, n, Float64[finalize_number(ctx, t) for t in e.seed]), 53
            steps = ["rounded to $T from the stored Float64 rule, without refinement; residual evaluated in Float64"]
        end
        guard, iterations = 0, 0
        cond = isnan(e.cond) ? 1.0 : e.cond                  # 1: not recorded
    else
        basis = moment_basis(e.structure, n)
        structure, θ64, seed_desc = seed_parameters(seed, e, ctx; family, domain, lower, basis, seed_desc)
        θ, res, guard = refine_seed(structure, n, θ64, ctx.bits; basis, cancel = ctx.cancel, verbose = ctx.verbose)
        res.converged || throw(RefinementError(family,
            "Gauss–Newton did not converge at degree $n (residual $(Float64(res.residual)) after $(res.iterations) " *
            "iterations); the seed is not returned unrefined"))
        m = margins(structure, θ)
        all(>(0), m) || throw(RefinementError(family,
            "the refined degree-$n rule left the admissible set (margins $(map(Float64, m)))"))
        wbits = ctx.bits + guard
        xs, ws = with_bits(() -> expand(structure, θ), wbits)
        # the certificate: the residual at the delivered parameters, at twice the precision
        θhat = [finalize_number(ctx, t) for t in θ]
        rbits = 2ctx.bits + guard
        resid = with_bits(() -> moment_residual(structure, n, BigFloat.(θhat), basis), rbits)
        steps = ["refine: Gauss–Newton with rank-revealing (column-pivoted QR) solve at $(wbits) bits",
                 @sprintf("guard: %d bits from measured cond(J) = %.2e", guard, res.cond_max)]
        cond, iterations = res.cond, res.iterations
    end
    D = dimension(domain)
    nodes = map(xs) do x
        p = cartesian_node(structure, x)
        SVector{D,T}(ntuple(j -> finalize_number(ctx, p[j]), D))
    end
    wt = [finalize_number(ctx, w) for w in ws]
    cert = Certificate(; equations, residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)), cond, iterations,
                       next_error = isfinite(e.next_error) ? BigFloat(e.next_error; precision = 64) : nothing)
    prov = Provenance(; family, derivation = Seeded(), path = vcat(["seed: " * seed_desc, "structure: " * string(structure)], steps),
                      seed_source = seed_desc, citations = vcat([CUBATURERULES_JL], citations), comparisons,
                      symmetry, license)
    return QuadratureRule(nodes, wt, domain, PolynomialDegree(n), prov, cert)
end
