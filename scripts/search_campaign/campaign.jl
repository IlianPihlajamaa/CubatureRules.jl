# One worker of the pooled search (see README.md):
#
#     julia --project --heap-size-hint=2G scripts/search_campaign/campaign.jl DOMAIN [options]
#
# DOMAIN is triangle, tetrahedron, square or cube. Options:
#     --root DIR        campaign directory, shared by all workers (default campaign/pooled)
#     --worker K        worker number, for file names, logs and the random stream (default 1)
#     --mode M          extend (default): work at one degree step above the best table, until
#                       --max-degree; then, and with improve, at a degree between --from and --to
#                       (default: the top four degrees of the best table)
#     --from N --to N   degree range for improve
#     --max-degree N
#     --blind N         ignore shipped entries at degree N and above (to time a rediscovery)
#     --hours H         wall-clock limit (default: none; a file named STOP in the root, or in
#                       the domain's directory, stops every worker in it)

using CubatureRules, LinearAlgebra, Printf, TOML, Dates, Random
const CR = CubatureRules
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "domains.jl"))
include(joinpath(@__DIR__, "engine.jl"))

function parse_options(args)
    opts = Dict{String,String}()
    i = 2
    while i <= length(args)
        startswith(args[i], "--") || error("unexpected argument $(args[i])")
        opts[args[i][3:end]] = args[i + 1]
        i += 2
    end
    return opts
end

if abspath(PROGRAM_FILE) == @__FILE__
    domain = ARGS[1]
    haskey(DOMAINS, domain) || error("domain must be one of $(join(sort(collect(keys(DOMAINS))), ", "))")
    o = parse_options(ARGS)
    getint(k) = haskey(o, k) ? parse(Int, o[k]) : nothing
    run_worker(domain, get(o, "root", joinpath(PKG, "campaign", "pooled")), something(getint("worker"), 1);
               mode = get(o, "mode", "extend"), lo = getint("from"), hi = getint("to"),
               maxdeg = something(getint("max-degree"), 1000),
               hours = haskey(o, "hours") ? parse(Float64, o["hours"]) : Inf,
               blind = something(getint("blind"), typemax(Int)))
end
