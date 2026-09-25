# A long-running search for more and smaller in-house rules. Meant to be left running.
#
#     julia --project -t 6 scripts/search_campaign.jl triangle|tetrahedron [maxdegree]
#
# Works on a copy of the shipped table in campaign/<domain>.toml, never on src/data: finds are
# merged into the shipped tables by hand, after scripts/certify_tables.jl has checked them.
# Progress goes to campaign/<domain>.log.
#
# It runs in rounds, each with fresh random seeds, forever:
#   1. extend: grow the highest degree in the table by one (src/refine/elimination.jl,
#      grow → eliminate chains), as far as maxdegree;
#   2. improve: for every degree already in the table, run new chains from the rule one
#      degree down, and keep a rule with fewer points than the entry has.
# Every rule recorded is refined to 100 digits first, and its seed is that refinement rounded
# to Float64, as the generators do.

using CubatureRules, LinearAlgebra, Printf, TOML, Dates
const CR = CubatureRules
BLAS.set_num_threads(1)

const DATA = joinpath(@__DIR__, "..", "src", "data")
const OUT = joinpath(@__DIR__, "..", "campaign")
const TABLES = Dict("triangle" => ("triangle_s3_seeds.toml", 3, 20),
                    "tetrahedron" => ("tetrahedron_s4_seeds.toml", 4, 8))
const CHAINS = 16

logline(io, args...) = (println(io, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "  ", args...); flush(io))

structure_of(e, N) = CR.SymmetricStructure([Int.(m) for m in e["structure"]], N)

function write_entries(path, header, entries)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => sort(entries; by = e -> e["degree"])); sorted = true)
    end
    mv(tmp, path; force = true)
end

"Grow the degree-(n−1) entry to degree n with `CHAINS` chains; a refined entry, or nothing."
function attempt(prev, n, N, seed, log)
    basis = CR.invariant_basis(N, n)
    t = @elapsed r = CR.grow_and_eliminate(structure_of(prev, N), Float64.(prev["seed"]), n;
                                           basis, chains = CHAINS, rng_seed = seed)
    if r === nothing
        logline(log, @sprintf("degree %2d: no chain succeeded (seed %#x, %.0f s)", n, seed, t))
        return nothing
    end
    s, θ0, counts = r
    θ, res, _ = CR.refine_symmetric(s, n, θ0, CR.digits_to_bits(100); basis)
    if !res.converged
        logline(log, @sprintf("degree %2d: %d points found but refinement failed", n, CR.npoints(s)))
        return nothing
    end
    θ64 = Float64.(θ)
    wmin, λmin = CR.rule_margins(s, θ64)
    logline(log, @sprintf("degree %2d: %d points (chains %s, seed %#x, %.0f s)", n, CR.npoints(s),
                          string(sort(counts)), seed, t))
    return Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => [o.mult for o in s.orbits],
                            "seed" => θ64, "status" => "ok", "note" => "", "method" => "campaign elimination from degree $(n - 1)",
                            "min_weight" => wmin, "min_barycentric" => λmin, "cond" => res.cond)
end

function main(domain, maxdeg)
    file, N, improve_from = TABLES[domain]
    mkpath(OUT)
    path = joinpath(OUT, file)
    isfile(path) || cp(joinpath(DATA, file), path)
    header = [l for l in eachline(path) if startswith(l, "#")]
    entries = TOML.parsefile(path)["rule"]
    log = open(joinpath(OUT, domain * ".log"), "a")
    logline(log, "campaign started: $domain up to degree $maxdeg, $(Threads.nthreads()) threads")
    round = 0
    while true
        round += 1
        byd = Dict(e["degree"] => e for e in entries if e["status"] == "ok")
        top = maximum(keys(byd))
        # 1. extend
        n = top + 1
        while n <= maxdeg
            e = attempt(byd[n - 1], n, N, UInt(0xca000000) + 1000 * round + n, log)
            e === nothing && break
            e["note"] = "found by scripts/search_campaign.jl"
            push!(entries, e)
            byd[n] = e
            write_entries(path, header, entries)
            logline(log, @sprintf("NEW degree %d: %d points", n, e["npoints"]))
            n += 1
        end
        # 2. improve, lowest degrees first (they are cheap, and a smaller rule there is a
        #    better starting point for the degrees above)
        for n in improve_from:maximum(keys(byd))
            haskey(byd, n - 1) || continue
            e = attempt(byd[n - 1], n, N, UInt(0xcb000000) + 1000 * round + n, log)
            e === nothing && continue
            if e["npoints"] < byd[n]["npoints"]
                old = byd[n]["npoints"]
                e["note"] = "found by scripts/search_campaign.jl; $(old - e["npoints"]) point(s) fewer than the previous entry"
                filter!(x -> x["degree"] != n, entries)
                push!(entries, e)
                byd[n] = e
                write_entries(path, header, entries)
                logline(log, @sprintf("IMPROVED degree %d: %d → %d points", n, old, e["npoints"]))
            end
        end
        logline(log, "round $round done")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 80)
end
