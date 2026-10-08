# Generate the Float64 seed tables for fully symmetric rules on the pyramid and the wedge,
# whose orbits are not those of the box: src/data/pyramid_c4v_seeds.toml and
# src/data/wedge_d3h_seeds.toml.
#
#     julia --project -t auto scripts/generate_orbit_seeds.jl pyramid 10 [chains]
#     julia --project -t auto scripts/generate_orbit_seeds.jl wedge 7:12 [chains]
#
# `n` searches degrees 1 to n, `a:b` degrees a to b, growing from the stored rule of degree
# a − 1. Stored entries outside the range are kept, and so is a stored entry with fewer points
# than the search finds.
#
# Every degree: neither domain is centrally symmetric. Up to `multistart` the point count is
# walked upwards from the previous degree's, with random starts on every orbit structure of
# each count whose unknowns exceed the equations by at most `slack` (refine/orbit_search.jl);
# above it each degree is grown from the one below and thinned by node elimination. No
# published numbers are used; the published point counts (Witherden & Vincent 2015) are only
# compared with in the log and the notes. Afterwards `scripts/certify_tables.jl pyramid wedge`
# refines every seed to 160 bits and records its residual.

using CubatureRules, Printf, TOML
const CR = CubatureRules

# name => (structure type, multistart up to this degree, slack, table, description, the
# structure and seed layout, published counts: Witherden & Vincent 2015, Table 1)
const DOMAINS = Dict(
    "pyramid" => (CR.PyramidStructure, 6, 4, "pyramid_c4v_seeds.toml",
                  "the pyramid with base [-1, 1]^2 at z = 0 and apex (0, 0, 1)",
                  ["`structure` lists the orbits by the sizes of the groups of equal nonzero collapsed",
                   "coordinates (ξ, η) = (x, y)/(1 - z); `seed` holds per orbit the weight of each point,",
                   "the distinct |ξ| values, then the height z."],
                  Dict(1 => 1, 2 => 5, 3 => 6, 4 => 10, 5 => 15, 6 => 24, 7 => 31, 8 => 47, 9 => 62, 10 => 83)),
    "wedge" => (CR.WedgeStructure, 6, 4, "wedge_d3h_seeds.toml",
                "the wedge, the triangle (0,0), (1,0), (0,1) times [-1, 1]",
                ["`structure` lists the orbits by the multiplicities of the triangle orbit's distinct",
                 "barycentric values, then 0 for an orbit at z = 0 or 2 for the pair ±z; `seed` holds per",
                 "orbit the weight of each point, the free barycentric values, then z for a pair."],
                Dict(1 => 1, 2 => 5, 3 => 8, 4 => 11, 5 => 16, 6 => 28, 7 => 35, 8 => 46, 9 => 60, 10 => 85)))

structure_toml(s::CR.PyramidStructure) = [o.base.mult for o in s.orbits]
structure_toml(s::CR.WedgeStructure) = CR.wedge_structure_lists(s)

function generate(name, degrees; chains = 16)
    T, multistart, slack, file, what, layout, published = DOMAINS[name]
    path = joinpath(@__DIR__, "..", "src", "data", file)
    old = isfile(path) ? TOML.parsefile(path)["rule"] : Dict{String,Any}[]
    stored = Dict(e["degree"] => e for e in old if e["status"] == "ok")
    # entries outside the range are kept as they are
    entries = filter(e -> !(e["degree"] in degrees), old)
    from(e) = (T([Int[x for x in m] for m in e["structure"]]), Float64.(e["seed"]))
    prev = haskey(stored, first(degrees) - 1) ? from(stored[first(degrees) - 1]) : nothing
    for n in degrees
        t = @elapsed got = if n <= multistart
            lo = prev === nothing ? 1 : CR.npoints(prev[1])
            CR.orbit_multistart(T, n; npts = lo:(4lo + 32), nstarts = 64, slack)
        elseif prev !== nothing
            r = CR.orbit_grow_and_eliminate(prev[1], prev[2], n; chains, rng_seed = 0xb0c5 + n)
            r === nothing ? nothing : (r[1], r[2])
        end
        # a stored rule with fewer points is kept
        if haskey(stored, n) && (got === nothing || stored[n]["npoints"] <= CR.npoints(got[1]))
            @printf("%s degree %2d: kept the stored %d-point rule%s [%.0f s]\n", name, n, stored[n]["npoints"],
                    got === nothing ? "" : " (found $(CR.npoints(got[1])))", t)
            push!(entries, stored[n])
            prev = from(stored[n])
            continue
        end
        if got === nothing
            @printf("%s degree %2d: nothing found [%.0f s]; stopping\n", name, n, t)
            append!(entries, (stored[k] for k in degrees if k >= n && haskey(stored, k)))
            break
        end
        s, θ = got
        w, c, _ = CR.structure_margins(s, θ)
        pub = get(published, n, nothing)
        @printf("%s degree %2d: %4d points (published %s), %2d orbits, min weight %.1e, distance from boundary %.3f [%.0f s]\n",
                name, n, CR.npoints(s), pub === nothing ? "—" : string(pub), length(s.orbits), w, c, t)
        flush(stdout)
        note = pub === nothing ? "beyond the degrees Witherden & Vincent (2015) tabulate" :
               CR.npoints(s) == pub ? "" :
               "$(abs(CR.npoints(s) - pub)) point(s) $(CR.npoints(s) > pub ? "more" : "fewer") than Witherden & Vincent (2015)"
        push!(entries, Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "structure" => structure_toml(s),
                                        "seed" => θ, "note" => note, "status" => "ok"))
        prev = (s, θ)
    end
    sort!(entries; by = e -> e["degree"])
    # the note scripts/certify_tables.jl adds to the header, for the entries kept
    header = isfile(path) ? collect(Iterators.takewhile(l -> startswith(l, "#"), eachline(path))) : String[]
    k = findfirst(l -> startswith(l, "# `residual`"), header)
    open(path, "w") do io
        println(io, "# Fully symmetric, positive-weight, interior rules on ", what,
                ", generated in-house by scripts/generate_orbit_seeds.jl.")
        foreach(l -> println(io, "# ", l), layout)
        k === nothing || foreach(l -> println(io, l), vcat("#", header[k:end]))
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    haskey(DOMAINS, ARGS[1]) || error("the domain is one of $(join(sort(collect(keys(DOMAINS))), ", "))")
    degrees = occursin(":", ARGS[2]) ? (:)(parse.(Int, split(ARGS[2], ":"))...) : 1:parse(Int, ARGS[2])
    generate(ARGS[1], degrees; chains = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16)
end
