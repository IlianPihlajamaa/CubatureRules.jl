# Generate the MIT-licensed seed table for fully symmetric tetrahedron rules,
# src/data/tetrahedron_s4_seeds.toml, from orbit structures alone (PLAN §0.2).
#
#     julia --project -t auto scripts/generate_tetrahedron_seeds.jl [mindegree] maxdegree
#
# No published point counts are relied on. Two searches are combined per degree, and the rule
# with the fewest points wins (ties: the most interior); an existing table entry is kept if
# neither search beats it.
#
#   elimination  16 chains of grow → eliminate starting from the previous degree's rule
#                (src/refine/elimination.jl). Cheap at every degree; gives an upper bound.
#   walk         multistart at increasing point counts, from a count below which no fully
#                symmetric rule exists up to the elimination bound, trying at each count the
#                structures closest in orbit-type mix to the previous degree's rule. Exhaustive
#                enough at low degree, but its cost grows 3–5× per degree, so it is only run
#                up to WALK_MAX_DEGREE.
#
# The counts are therefore the smallest these searches reached, not proven minima. The
# table is rewritten after every degree.

isdefined(Main, :CR) || (import CubatureRules; const CR = CubatureRules)
using LinearAlgebra, Printf, TOML
BLAS.set_num_threads(1)

const STRUCTURES_PER_COUNT = 6
const WALK_MAX_DEGREE = 10
const TYPES = sort(CR.orbit_pattern_types(4); by = p -> (-length(p.mult), p.mult))

"Number of orbits of each pattern type, in canonical type order."
type_counts(s) = [count(o -> o.mult == t.mult, s.orbits) for t in TYPES]
structure_of(e) = CR.SymmetricStructure([Int.(m) for m in e["structure"]], 4)

"Smallest point count with a structure having at least as many unknowns as equations."
function lower_count(n)
    npts = 1
    while isempty(CR.candidate_structures(npts, n, 4))
        npts += 1
    end
    return npts
end

function walk(n, counts, prevcounts, basis; nstarts)
    m = size(basis.Q, 2)
    for npts in counts
        ss = CR.candidate_structures(npts, n, 4; max_unknowns = m + 3)
        prevcounts === nothing ||
            sort!(ss; by = s -> (sum(abs, type_counts(s) .- prevcounts), CR.nunknowns(s)))  # stable
        for s in first(ss, STRUCTURES_PER_COUNT)
            t = @elapsed found = CR.multistart(s, n; nstarts, basis, first_only = false)
            @printf("    walk: degree %d, %d points, %d unknowns (%d equations): %d valid in %.0f s\n",
                    n, npts, CR.nunknowns(s), m, length(found), t)
            flush(stdout)
            isempty(found) && continue
            best = argmax(f -> (f[3], -f[4]), found)
            return s, best[1], "walk"
        end
    end
    return nothing
end

better(a, b) = b === nothing ||
               (CR.npoints(a[1]), -CR.rule_margins(a[1], a[2])[2]) < (CR.npoints(b[1]), -CR.rule_margins(b[1], b[2])[2])

function main(mindeg, maxdeg, path)
    old = isfile(path) ? TOML.parsefile(path)["rule"] : Dict{String,Any}[]
    entries = filter(e -> !(mindeg <= e["degree"] <= maxdeg), old)
    olds = Dict(e["degree"] => e for e in old if e["status"] == "ok")
    for n in mindeg:maxdeg
        basis = CR.invariant_basis(4, n)
        prev = get(olds, n - 1, nothing)
        best = nothing
        # 1. elimination from the previous degree's rule
        if prev !== nothing
            t = @elapsed r = CR.grow_and_eliminate(structure_of(prev), Float64.(prev["seed"]), n; basis, chains = 16)
            if r !== nothing
                best = (r[1], r[2], "elimination")
                @printf("    elimination: degree %d → %d points (chains %s) in %.0f s\n",
                        n, CR.npoints(r[1]), string(sort(r[3])), t)
                flush(stdout)
            end
        end
        # 2. the existing entry, if any
        if haskey(olds, n)
            cand = (structure_of(olds[n]), Float64.(olds[n]["seed"]), get(olds[n], "method", "walk"))
            better(cand, best) && (best = cand)
        end
        # 3. walk below the best count so far
        if n <= WALK_MAX_DEGREE
            lo = max(lower_count(n), prev === nothing ? 1 : prev["npoints"])
            hi = best === nothing ? lo + 40 : CR.npoints(best[1]) - 1
            w = walk(n, lo:hi, prev === nothing ? nothing : type_counts(structure_of(prev)), basis;
                     nstarts = n <= 8 ? 512 : 1024)
            w === nothing || !better(w, best) || (best = w)
        end
        if best === nothing
            @printf("degree %2d: no rule found\n", n)
            push!(entries, Dict("degree" => n, "npoints" => 0, "structure" => Vector{Int}[], "seed" => Float64[],
                                "status" => "failed", "note" => "no valid rule found"))
            write_table(path, entries)
            continue
        end
        s, θbest, method = best
        θ, res, _ = CR.refine_symmetric(s, n, θbest, CR.digits_to_bits(100); basis)
        res.converged || error("degree $n: refinement to 100 digits failed")
        θ64 = Float64.(θ)
        wmin, λmin = CR.rule_margins(s, θ64)
        @printf("degree %2d: %3d pts  %d unknowns  (%s)  cond %.1e  iters %d\n",
                n, CR.npoints(s), CR.nunknowns(s), method, res.cond, res.iterations)
        flush(stdout)
        e = Dict("degree" => n, "npoints" => CR.npoints(s), "structure" => [o.mult for o in s.orbits],
                 "seed" => θ64, "status" => "ok", "note" => "", "method" => method,
                 "min_weight" => wmin, "min_barycentric" => λmin, "cond" => res.cond)
        push!(entries, e)
        olds[n] = e
        write_table(path, entries)
    end
    return entries
end

function write_table(path, entries)
    open(path, "w") do io
        println(io, "# Fully symmetric (S₄) tetrahedron rule seeds for CubatureRules.jl — MIT licence.")
        println(io, "# Generated by scripts/generate_tetrahedron_seeds.jl from orbit structures alone;")
        println(io, "# no numbers or point counts were taken from any published table. See PROVENANCE.toml.")
        println(io, "#")
        println(io, "# seed layout: per orbit, [weight per point, free barycentric values...].")
        println(io)
        TOML.print(io, Dict("rule" => sort(entries; by = e -> e["degree"])); sorted = true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    lo, hi = length(ARGS) == 2 ? Tuple(parse.(Int, ARGS)) : (1, isempty(ARGS) ? 10 : parse(Int, ARGS[1]))
    main(lo, hi, joinpath(@__DIR__, "..", "src", "data", "tetrahedron_s4_seeds.toml"))
end
