# A long-running search for more and smaller fully symmetric rules on the square and the
# cube, on one core:
#
#     julia --project -t 1 scripts/search_box_campaign.jl <hours> [square_max cube_max]
#
# Works on copies of the shipped tables in campaign/, never on src/data: finds are merged into
# the shipped tables by hand, after scripts/certify_tables.jl has checked them. Each task is
# one grow → eliminate chain (src/refine/box_search.jl) from the rule two degrees down, with
# its own random seed. The tasks alternate between the square and the cube, and for each
# between extending the table by one odd degree (up to square_max, cube_max) and improving a
# degree already in it, keeping a rule with fewer points than the entry has. Every chain is
# logged in campaign/box_search.log, new and improved rules also in campaign/box.log.

using CubatureRules, LinearAlgebra, Printf, TOML, Dates, Random
const CR = CubatureRules
BLAS.set_num_threads(1)

const DATA = joinpath(@__DIR__, "..", "src", "data")
const OUT = joinpath(@__DIR__, "..", "campaign")
const FILES = Dict(2 => "square_d4_seeds.toml", 3 => "cube_oh_seeds.toml")
const IMPROVE_FROM = Dict(2 => 9, 3 => 7)

"A cancellation token that fires at a wall-clock deadline."
struct Deadline
    t::Float64
end
CR.iscancelled(d::Deadline) = time() > d.t

stamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
logline(path, msg) = open(io -> (println(io, stamp(), "  ", msg); flush(io)), path, "a")
read_entries(path) = [e for e in TOML.parsefile(path)["rule"]]
function write_entries(path, entries)
    header = [l for l in eachline(path) if startswith(l, "#")]
    tmp = path * ".tmp"
    open(tmp, "w") do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => sort(entries; by = e -> e["degree"])); sorted = true)
    end
    mv(tmp, path; force = true)
end
structure_of(e, D) = CR.BoxStructure([Int[x for x in m] for m in e["structure"]], D)

function main(hours, maxdeg)
    mkpath(OUT)
    for D in (2, 3)
        path = joinpath(OUT, FILES[D])
        isfile(path) || cp(joinpath(DATA, FILES[D]), path)
    end
    chainlog, mainlog = joinpath(OUT, "box_search.log"), joinpath(OUT, "box.log")
    deadline = time() + 3600hours
    rng = Random.Xoshiro(time_ns())
    logline(chainlog, "started for $hours h, up to degree $(maxdeg[2]) on the square and $(maxdeg[3]) on the cube")
    task = 0
    while time() < deadline
        task += 1
        D = isodd(task) ? 2 : 3
        path = joinpath(OUT, FILES[D])
        byd = Dict(e["degree"] => e for e in read_entries(path) if e["status"] == "ok")
        top = maximum(keys(byd))
        extend = iseven((task - 1) ÷ 2) && top < maxdeg[D]
        n = extend ? top + 2 : rand(rng, IMPROVE_FROM[D]:2:top)
        prev = get(byd, n - 2, nothing)
        prev === nothing && continue
        seed = rand(rng, UInt32)
        t = @elapsed got = try
            CR.box_grow_and_eliminate(structure_of(prev, D), Float64.(prev["seed"]), n; chains = 1, rng_seed = seed,
                                      cancel = Deadline(deadline))
        catch err
            err isa CR.CancelledError || rethrow()
            :stopped
        end
        got === :stopped && break
        name = D == 2 ? "square" : "cube"
        if got === nothing
            logline(chainlog, @sprintf("%-6s degree %2d: no rule (seed %#x, %.0f s)", name, n, seed, t))
            continue
        end
        s, θ, _ = got
        npts = CR.npoints(s)
        current = get(byd, n, nothing)
        logline(chainlog, @sprintf("%-6s degree %2d: %d points%s (seed %#x, %.0f s)", name, n, npts,
                                   current === nothing ? "" : " (table has $(current["npoints"]))", seed, t))
        current === nothing || npts < current["npoints"] || continue
        note = current === nothing ? "found by scripts/search_box_campaign.jl" :
               "found by scripts/search_box_campaign.jl; $(current["npoints"] - npts) point(s) fewer than the previous entry"
        entries = filter(e -> e["degree"] != n, read_entries(path))
        push!(entries, Dict{String,Any}("degree" => n, "npoints" => npts, "structure" => [o.mult for o in s.orbits],
                                        "seed" => θ, "status" => "ok", "note" => note))
        write_entries(path, entries)
        logline(mainlog, current === nothing ? "NEW $name degree $n: $npts points" :
                                               "IMPROVED $name degree $n: $(current["npoints"]) → $npts points")
    end
    logline(chainlog, "stopped after $task tasks")
end

if abspath(PROGRAM_FILE) == @__FILE__
    hours = parse(Float64, ARGS[1])
    main(hours, Dict(2 => length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 51, 3 => length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 29))
end
