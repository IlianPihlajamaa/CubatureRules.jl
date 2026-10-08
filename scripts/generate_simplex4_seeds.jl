# Generate the Float64 seed table for fully symmetric rules on the 4-simplex,
# src/data/simplex4_s5_seeds.toml, from orbit structures alone.
#
#     julia --project -t auto scripts/generate_simplex4_seeds.jl maxdegree [multistart] [chains]
#
# Up to degree `multistart` (default 4) the point count is walked upwards from the previous
# degree's, with random starts on every orbit structure of each count whose unknowns exceed
# the S₅-invariant equations by at most 3 (refine/seeds.jl); above it each degree is grown
# from the one below and thinned by node elimination (refine/elimination.jl), in `chains`
# independent chains. No published numbers are used. Afterwards
# `scripts/certify_tables.jl simplex4` refines every seed to 160 bits and records its residual.

using CubatureRules, Printf, TOML
const CR = CubatureRules

function generate(nmax; multistart = 4, chains = 16)
    entries = Dict{String,Any}[]
    prev = nothing
    for n in 1:nmax
        basis = CR.invariant_basis(5, n)
        m = size(basis.Q, 2)
        t = @elapsed got = if n <= multistart
            found = nothing
            lo = prev === nothing ? 1 : CR.npoints(prev[1])
            for npts in lo:(4lo + 64), s in CR.candidate_structures(npts, n, 5; max_unknowns = m + 3)
                f = CR.multistart(s, n; nstarts = 256, basis, first_only = false)
                isempty(f) && continue
                best = argmax(x -> (x[3], -x[4]), f)              # the most interior, then the earliest
                found = (s, best[1])
                break
            end
            found
        else
            r = CR.grow_and_eliminate(prev[1], prev[2], n; chains, basis, rng_seed = 0x5c + n)
            r === nothing ? nothing : (r[1], r[2])
        end
        if got === nothing
            @printf("4-simplex degree %2d: nothing found [%.0f s]; stopping\n", n, t)
            break
        end
        s, θ = got
        w, λ = CR.rule_margins(s, θ)
        @printf("4-simplex degree %2d: %4d points, %2d orbits, %d equations, min weight %.1e, min barycentric %.3f [%.0f s]\n",
                n, CR.npoints(s), length(s.orbits), m, w, λ, t)
        flush(stdout)
        push!(entries, Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s),
                                        "structure" => [o.mult for o in s.orbits], "seed" => θ, "status" => "ok"))
        prev = (s, θ)
    end
    path = joinpath(@__DIR__, "..", "src", "data", "simplex4_s5_seeds.toml")
    open(path, "w") do io
        println(io, "# Fully symmetric (S₅), positive-weight, interior rules on the 4-simplex, generated in-house")
        println(io, "# by scripts/generate_simplex4_seeds.jl. `structure` lists the orbits by the multiplicities of")
        println(io, "# their distinct barycentric coordinates; `seed` holds per orbit the weight of each point, then")
        println(io, "# all but the last distinct value.")
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate(parse(Int, ARGS[1]); multistart = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4,
             chains = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16)
end
