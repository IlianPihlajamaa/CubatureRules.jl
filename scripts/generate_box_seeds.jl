# Generate the Float64 seed tables for fully symmetric rules on the square, the cube and the
# 4-cube, and on the disk and the ball: src/data/square_d4_seeds.toml, cube_oh_seeds.toml,
# hypercube_b4_seeds.toml, disk_d4_seeds.toml and ball_oh_seeds.toml.
#
#     julia --project -t auto scripts/generate_box_seeds.jl square 31
#     julia --project -t auto scripts/generate_box_seeds.jl cube 15 [chains]
#     julia --project -t auto scripts/generate_box_seeds.jl hypercube 9 [chains]
#     julia --project -t auto scripts/generate_box_seeds.jl disk 21 [chains]
#     julia --project -t auto scripts/generate_box_seeds.jl ball 15 [chains]
#
# The orbits are those of the signed permutations of the coordinates in each case
# (symmetry/box.jl); the box and the ball differ in their moment systems and in where an
# orbit's values may lie (refine/box_search.jl).
#
# Odd degrees only: every orbit is centrally symmetric, so a rule of degree 2k is one of
# degree 2k + 1. Up to `multistart` the point count is walked upwards from the previous
# degree's, with random starts on every orbit structure of each count; above it each degree is
# grown from the one below and thinned by node elimination. No published numbers are used;
# the published point counts (Witherden & Vincent 2015) are only compared with in the log.
# Afterwards `scripts/certify_tables.jl square cube` refines every seed to 160 bits and
# records its residual.

using CubatureRules, Printf, TOML
const CR = CubatureRules

# name => (dimension, geometry, multistart up to this degree, table, description, published counts)
const DOMAINS = Dict(
    "square" => (2, CR.BoxGeometry(), 13, "square_d4_seeds.toml", "the square [-1, 1]^2",
                 Dict(1 => 1, 3 => 4, 5 => 8, 7 => 12, 9 => 20, 11 => 28, 13 => 37, 15 => 48, 17 => 60, 19 => 72, 21 => 85)),
    "cube" => (3, CR.BoxGeometry(), 7, "cube_oh_seeds.toml", "the cube [-1, 1]^3",
               Dict(1 => 1, 3 => 6, 5 => 14, 7 => 34, 9 => 58, 11 => 90)),
    "hypercube" => (4, CR.BoxGeometry(), 5, "hypercube_b4_seeds.toml", "the 4-cube [-1, 1]^4", Dict{Int,Int}()),
    "disk" => (2, CR.RoundGeometry(), 13, "disk_d4_seeds.toml", "the unit disk", Dict{Int,Int}()),
    "ball" => (3, CR.RoundGeometry(), 7, "ball_oh_seeds.toml", "the unit ball", Dict{Int,Int}()))

function generate(name, nmax; chains = 16)
    D, geom, multistart, file, what, published = DOMAINS[name]
    entries = Dict{String,Any}[]
    prev = nothing
    for n in 1:2:nmax
        t = @elapsed got = if n <= multistart
            lo = prev === nothing ? 1 : CR.npoints(prev[1])
            CR.box_multistart(n, D; npts = lo:(4lo + 64), nstarts = 64, geom)
        else
            r = CR.box_grow_and_eliminate(prev[1], prev[2], n; chains, rng_seed = 0xb0c5 + n, geom)
            r === nothing ? nothing : (r[1], r[2])
        end
        if got === nothing
            @printf("%s degree %2d: nothing found [%.0f s]; stopping\n", name, n, t)
            break
        end
        s, θ = got
        w, c, g = CR.orbit_margins(geom, s, θ)
        pub = get(published, n, nothing)
        @printf("%s degree %2d: %4d points (published %s), %2d orbits, min weight %.1e, distance from boundary %.3f [%.0f s]\n",
                name, n, CR.npoints(s), pub === nothing ? "—" : string(pub), length(s.orbits), w, c, t)
        flush(stdout)
        note = isempty(published) ? "" : pub === nothing ? "beyond the degrees Witherden & Vincent (2015) tabulate" :
               CR.npoints(s) == pub ? "" :
               "$(abs(CR.npoints(s) - pub)) point(s) $(CR.npoints(s) > pub ? "more" : "fewer") than Witherden & Vincent (2015)"
        push!(entries, Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s),
                                        "structure" => [o.mult for o in s.orbits], "seed" => θ,
                                        "note" => note, "status" => "ok"))
        prev = (s, θ)
    end
    path = joinpath(@__DIR__, "..", "src", "data", file)
    open(path, "w") do io
        println(io, "# Fully symmetric, positive-weight, interior rules on ", what,
                ", generated in-house by scripts/generate_box_seeds.jl.")
        println(io, "# `structure` lists the orbits by the sizes of their groups of equal nonzero coordinates;")
        println(io, "# `seed` holds per orbit the weight of each point, then the distinct coordinate values.")
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    haskey(DOMAINS, ARGS[1]) || error("the domain is one of $(join(sort(collect(keys(DOMAINS))), ", "))")
    generate(ARGS[1], parse(Int, ARGS[2]); chains = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16)
end
