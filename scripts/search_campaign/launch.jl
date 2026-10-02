# Start pooled-search workers as independent processes (see README.md):
#
#     julia --project scripts/search_campaign/launch.jl [--root DIR] SPEC [SPEC ...]
#
# Each SPEC is DOMAIN:MODE:COUNT[:FROM-TO], e.g.
#
#     triangle:extend:2 triangle:improve:1 tetrahedron:extend:2 square:extend:1 cube:extend:1
#
# Workers are single-threaded, run until a STOP file appears in the root, and keep running
# when this launcher exits. Their process ids are appended to <root>/workers.txt.

const ROOT_DEFAULT = normpath(joinpath(@__DIR__, "..", "..", "campaign", "pooled"))

function launch(root, specs)
    mkpath(root)
    julia = Base.julia_cmd()
    project = normpath(joinpath(@__DIR__, "..", ".."))
    script = joinpath(@__DIR__, "campaign.jl")
    next_id = Dict{String,Int}()
    for spec in specs
        parts = split(spec, ":")
        domain, mode, count = parts[1], parts[2], parse(Int, parts[3])
        range = length(parts) >= 4 ? split(parts[4], "-") : String[]
        for _ in 1:count
            id = (next_id[domain] = get(next_id, domain, existing_workers(root, domain)) + 1)
            args = ["--root", root, "--worker", string(id), "--mode", mode]
            isempty(range) || append!(args, ["--from", range[1], "--to", range[2]])
            mkpath(joinpath(root, domain, "logs"))
            out = joinpath(root, domain, "logs", "worker$(id).out")
            cmd = pipeline(`$julia --project=$project --threads=1 --heap-size-hint=2G $script $domain $args`;
                           stdout = out, stderr = out)
            p = run(detach(cmd); wait = false)
            pid = getpid(p)
            open(io -> println(io, "$domain $mode worker $id pid $pid"), joinpath(root, "workers.txt"), "a")
            println("started $domain $mode worker $id (pid $pid)")
        end
    end
end

"Worker numbers already used for a domain, so a second launch does not reuse them."
function existing_workers(root, domain)
    dir = joinpath(root, domain, "logs")
    isdir(dir) || return 0
    ids = [parse(Int, m[1]) for f in readdir(dir) for m in (match(r"^worker(\d+)\.log$", f),) if m !== nothing]
    return isempty(ids) ? 0 : maximum(ids)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = copy(ARGS)
    root = ROOT_DEFAULT
    if !isempty(args) && args[1] == "--root"
        root = args[2]; args = args[3:end]
    end
    launch(root, args)
end
