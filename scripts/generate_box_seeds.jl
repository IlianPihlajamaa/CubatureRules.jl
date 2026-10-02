# Generate the Float64 seed tables for fully symmetric rules on the square and the cube,
# src/data/square_d4_seeds.toml and src/data/cube_oh_seeds.toml.
#
#     julia --project -t auto scripts/generate_box_seeds.jl square 31
#     julia --project -t auto scripts/generate_box_seeds.jl cube 15 [chains]
#
# Odd degrees only: every orbit is centrally symmetric, so a rule of degree 2k is one of
# degree 2k + 1. Up to `MULTISTART_MAX` the point count is walked upwards from the previous
# degree's, with random starts on every orbit structure of each count; above it each degree is
# grown from the one below and thinned by node elimination. No published numbers are used;
# the published point counts (Witherden & Vincent 2015) are only compared with in the log.
# Afterwards `scripts/certify_tables.jl square cube` refines every seed to 160 bits and
# records its residual.

using CubatureRules, Printf, TOML
const CR = CubatureRules

const MULTISTART_MAX = Dict(2 => 13, 3 => 7)
const PUBLISHED = Dict(2 => Dict(1 => 1, 3 => 4, 5 => 8, 7 => 12, 9 => 20, 11 => 28, 13 => 37, 15 => 48, 17 => 60,
                                 19 => 72, 21 => 85),
                       3 => Dict(1 => 1, 3 => 6, 5 => 14, 7 => 34, 9 => 58, 11 => 90))
const FILES = Dict(2 => "square_d4_seeds.toml", 3 => "cube_oh_seeds.toml")

function generate(D, nmax; chains = 16)
    entries = Dict{String,Any}[]
    prev = nothing
    for n in 1:2:nmax
        t = @elapsed got = if n <= MULTISTART_MAX[D]
            lo = prev === nothing ? 1 : CR.npoints(prev[1])
            CR.box_multistart(n, D; npts = lo:(4lo + 64), nstarts = 64)
        else
            r = CR.box_grow_and_eliminate(prev[1], prev[2], n; chains, rng_seed = 0xb0c5 + n)
            r === nothing ? nothing : (r[1], r[2])
        end
        if got === nothing
            @printf("D=%d degree %2d: nothing found [%.0f s]; stopping\n", D, n, t)
            break
        end
        s, θ = got
        w, c, g = CR.box_margins(s, θ)
        pub = get(PUBLISHED[D], n, nothing)
        @printf("D=%d degree %2d: %4d points (published %s), %2d orbits, min weight %.1e, distance from boundary %.3f [%.0f s]\n",
                D, n, CR.npoints(s), pub === nothing ? "—" : string(pub), length(s.orbits), w, c, t)
        flush(stdout)
        note = pub === nothing ? "beyond the degrees Witherden & Vincent (2015) tabulate" :
               CR.npoints(s) == pub ? "" :
               "$(abs(CR.npoints(s) - pub)) point(s) $(CR.npoints(s) > pub ? "more" : "fewer") than Witherden & Vincent (2015)"
        push!(entries, Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s),
                                        "structure" => [o.mult for o in s.orbits], "seed" => θ,
                                        "note" => note, "status" => "ok"))
        prev = (s, θ)
    end
    path = joinpath(@__DIR__, "..", "src", "data", FILES[D])
    open(path, "w") do io
        println(io, "# Fully symmetric, positive-weight, interior rules on the ", D == 2 ? "square" : "cube",
                " [-1, 1]^", D, ", generated in-house by scripts/generate_box_seeds.jl.")
        println(io, "# `structure` lists the orbits by the sizes of their groups of equal nonzero coordinates;")
        println(io, "# `seed` holds per orbit the weight of each point, then the distinct coordinate values.")
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    D = ARGS[1] == "square" ? 2 : ARGS[1] == "cube" ? 3 : error("square or cube")
    generate(D, parse(Int, ARGS[2]); chains = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16)
end
