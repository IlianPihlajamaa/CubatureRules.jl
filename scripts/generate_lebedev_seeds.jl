# Generate the O_h-symmetric sphere seed table in src/data/lebedev_seeds.toml.
#
#     julia --project -t auto scripts/generate_lebedev_seeds.jl [mindeg] [maxdeg]
#
# For each odd degree it walks point counts upward and, at each count, tries every orbit
# structure with at least as many unknowns as invariant equations, keeping the first that
# yields a rule that is positive, non-degenerate and has distinct nodes. Nothing is copied
# from any published table: only the counts are afterwards compared, in the log.
#
# Even degrees are skipped on purpose. Every O_h orbit is centrally symmetric, so a rule of
# degree 2k is automatically of degree 2k+1: the odd harmonics integrate to zero whatever
# the parameters.

using CubatureRules, Printf, TOML, Dates
const CR = CubatureRules

# Lebedev's own point counts, for comparison in the log only — never used as input.
const LEBEDEV_COUNTS = Dict(3 => 6, 5 => 14, 7 => 26, 9 => 38, 11 => 50, 13 => 74, 15 => 86,
                            17 => 110, 19 => 146, 21 => 170, 23 => 194, 25 => 230, 27 => 266,
                            29 => 302, 31 => 350, 35 => 434, 41 => 590, 47 => 770)

const SEED_FILE = joinpath(@__DIR__, "..", "src", "data", "lebedev_seeds.toml")

"Try one structure; return (θ, wmin, dmin) or nothing."
function try_structure(s, n, nstarts, rng_seed)
    found = CR.multistart_octahedral(s, n; nstarts, rng_seed, first_only = false)
    isempty(found) && return nothing
    best = argmax(t -> (t[2], t[3]), found)        # largest minimum weight, then margin
    return best[1], best[2], best[3]
end

function search_degree(n; maxpts, nstarts)
    m = length(CR.invariant_exponents(n))
    for npts in 6:2:maxpts
        for s in CR.octahedral_candidate_structures(npts, n; max_unknowns = m + 4)
            got = try_structure(s, n, nstarts, 0x5eed + n)
            got === nothing && continue
            return (npts, s, got...)
        end
    end
    return nothing
end

function write_table(path, entries)
    open(path, "w") do io
        println(io, "# O_h-symmetric sphere rules: orbit structures and Float64 seeds.")
        println(io, "# Generated in-house by scripts/generate_lebedev_seeds.jl; no published")
        println(io, "# numbers were used. See src/data/PROVENANCE.toml.")
        println(io, "# Last written: ", Dates.format(Dates.now(), "yyyy-mm-dd"))
        for e in sort(entries; by = x -> x["degree"])
            println(io, "\n[[rule]]")
            for k in sort(collect(keys(e)))
                v = e[k]
                print(io, k, " = ")
                if v isa String
                    println(io, repr(v))
                elseif v isa AbstractVector{<:AbstractString}
                    println(io, "[", join(repr.(v), ", "), "]")
                elseif v isa AbstractVector
                    println(io, "[", join(string.(v), ", "), "]")
                else
                    println(io, v)
                end
            end
        end
    end
end

function main(mindeg, maxdeg; maxpts = 600, nstarts = 2048)
    entries = isfile(SEED_FILE) ?
              filter(e -> !(mindeg <= e["degree"] <= maxdeg), TOML.parsefile(SEED_FILE)["rule"]) :
              Dict{String,Any}[]
    for n in mindeg:2:maxdeg
        isodd(n) || continue
        t0 = time()
        got = search_degree(n; maxpts, nstarts)
        if got === nothing
            @printf("degree %2d: nothing found up to %d points (%.0f s)\n", n, maxpts, time() - t0)
            flush(stdout)
            continue
        end
        npts, s, θ, wmin, dmin = got
        θ, res, guard = CR.refine_octahedral(s, n, θ, CR.digits_to_bits(100))
        res.converged || error("degree $n: refinement to 100 digits failed")
        θ64 = Float64.(θ)
        published = get(LEBEDEV_COUNTS, n, nothing)
        note = published === nothing ? "" :
               npts == published ? "matches the published Lebedev count" :
               npts < published ? "$(published - npts) point(s) fewer than the published count" :
               "$(npts - published) point(s) more than the published count, which is not positive"
        @printf("degree %2d: %3d points %-30s (Lebedev %s)  wmin %.2e  margin %.2e  cond %.1e  (%.0f s)\n",
                n, npts, string([o.kind for o in s.orbits]),
                published === nothing ? "—" : string(published), wmin, dmin, res.cond, time() - t0)
        flush(stdout)
        push!(entries, Dict{String,Any}(
            "degree" => n, "npoints" => npts, "structure" => [string(o.kind) for o in s.orbits],
            "seed" => θ64, "status" => "ok", "note" => note, "method" => "multistart",
            "min_weight" => wmin, "min_margin" => dmin, "cond" => res.cond))
        write_table(SEED_FILE, entries)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mindeg = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3
    maxdeg = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 29
    main(mindeg, maxdeg)
end
