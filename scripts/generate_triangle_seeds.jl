# Generate the MIT-licensed seed table for fully symmetric triangle rules,
# src/data/triangle_s3_seeds.toml, from orbit structures alone (PLAN §0.2).
#
#     julia --project -t auto scripts/generate_triangle_seeds.jl [mindegree] maxdegree
#
# Entries for degrees in the range are replaced; the others are kept. The table is rewritten
# after every degree, so a long run can be interrupted without losing finished degrees.
#
# For each degree the target point count is the minimal count reported by Xiao & Gimbutas
# (2010) — a combinatorial fact stated in the paper, not tabulated data. Orbit structures
# with that many points and at least as many unknowns as S₃-invariant equations are searched
# by multistart Levenberg–Marquardt; the first structure (fewest unknowns) with a valid rule
# wins, and among its valid rules the one with the largest minimum barycentric coordinate
# is kept. The Float64 seed written out is the 100-digit refinement rounded to Float64. If
# no structure at the target count yields a rule, the count is increased and the entry
# says so.
#
# Above degree 20, where multistart can afford only a few structures per count, node
# elimination from the previous degree's rule runs first (src/refine/elimination.jl), and
# multistart is only tried at counts below the one elimination reached.

isdefined(Main, :CR) || (import CubatureRules; const CR = CubatureRules)
using LinearAlgebra, Printf, TOML
BLAS.set_num_threads(1)

# Minimal point counts of fully symmetric rules, Xiao & Gimbutas (2010), Table 1, column n6.
const XG_POINTS = [1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79,
                   87, 96, 103, 112, 120, 130, 141, 150, 159, 171, 181, 193, 204, 214, 228,
                   243, 252, 267, 282, 295, 309, 324, 339, 354, 370, 385, 399, 423, 435, 453]

const ELIM_FROM_DEGREE = 20

function search(n, npts; nstarts, maxextra = 3)
    basis = CR.invariant_basis(3, n)
    # structures with the fewest unknowns first; at high degree only a few, since each one
    # costs thousands of starts
    maxstruct = n <= 20 ? typemax(Int) : 3
    for extra in 0:maxextra, s in first(CR.candidate_structures(npts + extra, n), maxstruct)
        t = @elapsed found = CR.multistart(s, n; nstarts, basis, first_only = false)
        @printf("    degree %d, %d points, %d unknowns: %d valid in %.0f s\n",
                n, npts + extra, CR.nunknowns(s), length(found), t)
        flush(stdout)
        isempty(found) && continue
        best = argmax(t -> (t[3], -t[4]), found)
        return s, best[1], "multistart ($(length(found)) valid)"
    end
    return nothing
end

"Grow the degree-(n-1) entry and eliminate points (src/refine/elimination.jl)."
function eliminate_from(prev, n)
    s = CR.SymmetricStructure([Int.(m) for m in prev["structure"]], 3)
    t = @elapsed r = CR.grow_and_eliminate(s, Float64.(prev["seed"]), n; chains = 16)
    r === nothing && return nothing
    @printf("    elimination: degree %d → %d points (chains %s) in %.0f s\n", n, CR.npoints(r[1]), string(sort(r[3])), t)
    flush(stdout)
    return r[1], r[2], "elimination from degree $(n - 1)"
end

function main(mindeg, maxdeg, path)
    entries = isfile(path) ? filter(e -> !(mindeg <= e["degree"] <= maxdeg), TOML.parsefile(path)["rule"]) :
              Dict{String,Any}[]
    for n in mindeg:maxdeg
        npts = XG_POINTS[n]
        nstarts = n <= 10 ? 512 : n <= 15 ? 1024 : n <= 20 ? 2048 : 4096
        prev = findfirst(e -> e["degree"] == n - 1 && e["status"] == "ok", entries)
        # above degree 20 the multistart search only affords a few structures per count, so
        # node elimination runs first; multistart is then only tried below its count
        elim = n > ELIM_FROM_DEGREE && prev !== nothing ? eliminate_from(entries[prev], n) : nothing
        best = elim
        if elim === nothing || CR.npoints(elim[1]) > npts
            maxextra = elim === nothing ? 3 : CR.npoints(elim[1]) - npts - 1
            found = search(n, npts; nstarts, maxextra)
            found === nothing || (best = found)
        end
        if best === nothing
            @printf("degree %2d: no rule found near %d points\n", n, npts)
            push!(entries, Dict("degree" => n, "npoints" => npts, "structure" => Vector{Int}[],
                                "seed" => Float64[], "status" => "failed",
                                "note" => "no valid rule found within 3 points of the published count"))
            write_table(path, entries)
            continue
        end
        s, θ0, method = best
        θ, res, guard = CR.refine_symmetric(s, n, θ0, CR.digits_to_bits(100); basis = CR.invariant_basis(3, n))
        res.converged || error("degree $n: refinement to 100 digits failed")
        θ64 = Float64.(θ)
        wmin, λmin = CR.rule_margins(s, θ64)
        extra = CR.npoints(s) - npts
        @printf("degree %2d: %3d pts (published %d)  %s  cond %.1e  iters %d\n",
                n, CR.npoints(s), npts, method, res.cond, res.iterations)
        flush(stdout)
        note = extra == 0 ? "" : extra > 0 ? "no valid rule found at the published count; $(extra) extra point(s)" :
               "$(-extra) point(s) fewer than the published count"
        push!(entries, Dict("degree" => n, "npoints" => CR.npoints(s),
                            "structure" => [o.mult for o in s.orbits], "seed" => θ64,
                            "status" => "ok", "note" => note, "method" => method,
                            "min_weight" => wmin, "min_barycentric" => λmin, "cond" => res.cond))
        write_table(path, entries)
    end
    return entries
end

function write_table(path, entries)
    open(path, "w") do io
        println(io, "# Fully symmetric (S₃) triangle rule seeds for CubatureRules.jl — MIT licence.")
        println(io, "# Generated by scripts/generate_triangle_seeds.jl from orbit structures alone;")
        println(io, "# no numbers were taken from any published table. See PROVENANCE.toml.")
        println(io, "#")
        println(io, "# seed layout: per orbit, [weight per point, free barycentric values...].")
        println(io)
        TOML.print(io, Dict("rule" => sort(entries; by = e -> e["degree"])); sorted = true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    lo, hi = length(ARGS) == 2 ? Tuple(parse.(Int, ARGS)) : (1, isempty(ARGS) ? 20 : parse(Int, ARGS[1]))
    main(lo, hi, joinpath(@__DIR__, "..", "src", "data", "triangle_s3_seeds.toml"))
end
