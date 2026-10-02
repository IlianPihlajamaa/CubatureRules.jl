# Summarise a pooled search (see README.md):
#
#     julia --project scripts/search_campaign/status.jl [ROOT]
#
# For every domain under ROOT: the degrees where the campaign has beaten or extended the shipped
# table, the pool at the top degrees (parents by distance above square), and recent finds.

using TOML, Printf

const PKG = normpath(joinpath(@__DIR__, "..", ".."))
const FILES = Dict("triangle" => "triangle_s3_seeds.toml", "tetrahedron" => "tetrahedron_s4_seeds.toml",
                   "square" => "square_d4_seeds.toml", "cube" => "cube_oh_seeds.toml")

function status(root)
    for domain in sort(collect(keys(FILES)))
        dir = joinpath(root, domain)
        isdir(dir) || continue
        shipped = Dict(e["degree"] => e["npoints"] for e in TOML.parsefile(joinpath(PKG, "src", "data", FILES[domain]))["rule"])
        best = Dict{Int,Int}()
        exc = Dict{Int,Int}()
        cdir = joinpath(dir, "candidates")
        if isdir(cdir)
            for f in readdir(cdir; join = true)
                endswith(f, ".toml") || continue
                e = TOML.parsefile(f)
                n = e["degree"]
                e["npoints"] < get(best, n, typemax(Int)) || continue
                best[n] = e["npoints"]; exc[n] = get(e, "excess", 0)
            end
        end
        println("== ", domain, ": shipped up to degree ", maximum(keys(shipped)))
        for n in sort(collect(keys(best)))
            s = get(shipped, n, nothing)
            println(@sprintf("   degree %3d: %5d points", n, best[n]),
                    exc[n] > 0 ? "  excess $(exc[n])" : "",
                    s === nothing ? "  (new)" : s > best[n] ? "  (shipped $s)" : "  (shipped $s, not better)")
        end
        pdir = joinpath(dir, "pool")
        if isdir(pdir)
            for d in sort(readdir(pdir))[max(1, end - 2):end]
                names = filter(f -> endswith(f, ".toml"), readdir(joinpath(pdir, d)))
                ex = [parse(Int, m[1]) for f in names for m in (match(r"^e(\d+)_", f),) if m !== nothing]
                hist = join(("$(e):$(count(==(e), ex))" for e in sort(unique(ex))), " ")
                println("   pool ", d, ": ", length(names), " parents (excess:count ", hist, ")")
            end
        end
        logs = joinpath(dir, "logs")
        if isdir(logs)
            events = String[]
            for f in readdir(logs; join = true)
                endswith(f, ".log") || continue
                append!(events, filter(l -> occursin(r"NEW degree|IMPROVED|error", l), readlines(f)))
            end
            foreach(l -> println("   ", l), sort(events)[max(1, end - 4):end])
        end
    end
end

status(isempty(ARGS) ? joinpath(PKG, "campaign", "pooled") : ARGS[1])
