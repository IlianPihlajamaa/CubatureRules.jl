# The pooled search engine (see README.md): grow / prune / finish workers around a shared pool
# of exact rules on disk. Domain-specific work goes through the adapters in domains.jl.

const PKG = normpath(joinpath(@__DIR__, "..", ".."))
const FINISH_EXCESS = 3                    # parents at most this far above square are "near square"
grow_excess(m, rng) = rand(rng, (3, 6, 12, max(8, m ÷ 8), max(8, m ÷ 4)))   # how far above square a grown rule starts
const SIGMAS = (0.0, 1e-4, 5e-4, 2e-3, 6e-3, 1.5e-2)   # coordinate perturbations in the finish step
const CHAIN_SHARE = 0.25                   # share of a worker's time spent in the package's own chains
const PRUNE_TRIES = 10_000                 # refits per prune step: in effect all, as `eliminate` does
const MIN_SEPARATION = 1e-6                # as scripts/verify_tables.jl: closer nodes are merging orbits
const PATIENCE = 2.0                       # hours a worker spends on an unsettled frontier rule before moving up

mutable struct Worker
    d::CampaignDomain
    root::String                           # <campaign root>/<domain>
    id::Int
    rng::MersenneTwister
    log::IO
    counter::Int
    tries::Dict{String,Int}                # finish attempts per parent file, by this worker
    family_tries::Dict{Any,Int}
    families::Dict{String,Any}             # parent file → family key
    best::Dict{Int,Int}                    # degree → fewest points known
    excess::Dict{Int,Int}                  # degree → excess of that rule (0 for shipped and imported rules)
    worked::Dict{Tuple{Int,Int},Float64}   # (degree, best points) → seconds this worker spent there since
    shipped_top::Int
    best_read::Float64
    stats::Dict{Symbol,Int}
    blind::Int                             # shipped entries at this degree and above are ignored (benchmarks)
    exhausted::Set{String}                 # pool files no prune move improved: left to the finish step
end

stamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
say(w::Worker, msg) = (println(w.log, stamp(), "  ", msg); flush(w.log))

pooldir(w, n) = joinpath(w.root, "pool", @sprintf("deg%03d", n))
candir(w) = joinpath(w.root, "candidates")

"Write `path` through a temporary file and a rename, so readers never see a partial file."
function atomic_write(f, path)
    mkpath(dirname(path))
    tmp = path * ".tmp" * string(getpid())
    open(f, tmp, "w")
    mv(tmp, path; force = true)
end

function save_rule(w, path, n, s, θ, extra = Dict{String,Any}())
    d = Dict{String,Any}("degree" => n, "npoints" => CR.npoints(s), "unknowns" => CR.nunknowns(s),
                         "equations" => n_equations(w.d, n), "structure" => orbit_mults(w.d, s),
                         "seed" => collect(Float64, θ), "worker" => w.id, "found" => stamp())
    merge!(d, extra)
    atomic_write(io -> TOML.print(io, d; sorted = true), path)
end
function load_rule(w, path)
    t = TOML.parsefile(path)
    return t["degree"], structure(w.d, t["structure"]), Float64.(t["seed"]), t
end

const POOLNAME = r"^e(\d+)_p(\d+)_h([0-9a-f]+)\.toml$"
function pool_path(w, n, s, θ)
    h = string(hash((orbit_mults(w.d, s), round.(θ; digits = 7))); base = 16)
    return joinpath(pooldir(w, n), @sprintf("e%03d_p%06d_h%s.toml", CR.nunknowns(s) - n_equations(w.d, n), CR.npoints(s), h))
end
"The pool at degree n: (path, excess, points), from the file names alone."
function pool_list(w, n)
    dir = pooldir(w, n)
    isdir(dir) || return NamedTuple{(:path, :excess, :npoints),Tuple{String,Int,Int}}[]
    out = NamedTuple{(:path, :excess, :npoints),Tuple{String,Int,Int}}[]
    for f in readdir(dir)
        m = match(POOLNAME, f)
        m === nothing || push!(out, (path = joinpath(dir, f), excess = parse(Int, m[1]), npoints = parse(Int, m[2])))
    end
    return out
end

# --- the best table: shipped entries plus every candidate found --------------------------------

shipped(w) = [e for e in TOML.parsefile(joinpath(PKG, "src", "data", table_file(w.d)))["rule"]
              if get(e, "status", "ok") == "ok" && e["degree"] < w.blind]
function best_entries(w)
    best = Dict(e["degree"] => e for e in shipped(w))
    isdir(candir(w)) || return best
    for f in readdir(candir(w); join = true)
        endswith(f, ".toml") || continue
        e = try TOML.parsefile(f) catch; continue end
        n = e["degree"]
        (!haskey(best, n) || e["npoints"] < best[n]["npoints"]) && (best[n] = e)
    end
    return best
end
function refresh_best!(w; force = false)
    (force || time() - w.best_read > 60) || return w.best
    es = best_entries(w)
    w.best = Dict(n => e["npoints"] for (n, e) in es)
    w.excess = Dict(n => get(e, "excess", 0) for (n, e) in es)
    w.best_read = time()
    return w.best
end

"Write the merged best table, in the shipped format, as <root>/<domain>/<table file>."
function write_best_table(w)
    header = [l for l in eachline(joinpath(PKG, "src", "data", table_file(w.d))) if startswith(l, "#")]
    entries = [Dict(k => v for (k, v) in e if !(k in ("found", "worker", "parent", "how", "excess"))) for e in values(best_entries(w))]
    atomic_write(joinpath(w.root, table_file(w.d))) do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => sort(entries; by = e -> e["degree"])); sorted = true)
    end
end

"Smallest max-norm distance between two nodes, as scripts/verify_tables.jl measures it."
function min_separation(s, θ)
    x = try
        first(CR.expand(s, Float64.(θ)))
    catch                                     # a box orbit whose points coincide
        return 0.0
    end
    sep = Inf
    @inbounds for i in eachindex(x), j in (i + 1):length(x)
        dmax = 0.0
        for k in eachindex(x[i])
            dmax = max(dmax, abs(x[i][k] - x[j][k]))
        end
        sep = min(sep, dmax)
    end
    return sep
end

# --- recording what a solve found --------------------------------------------------------------

"""
Keep a valid degree-`n` rule: in the pool if it still has unknowns to spare, and as a candidate,
refined to 100 digits, if it has fewer points than any rule known at that degree.
"""
function record!(w, n, s, θ, how)
    excess = CR.nunknowns(s) - n_equations(w.d, n)
    npts = CR.npoints(s)
    if excess >= 1
        # one file per rule: the name carries a hash of the structure and the rounded parameters,
        # so a rule found again (often, by several workers) is not stored twice
        path = pool_path(w, n, s, θ)
        isfile(path) || save_rule(w, path, n, s, θ, Dict{String,Any}("how" => how))
    end
    npts < get(refresh_best!(w), n, typemax(Int)) || return            # the cached counts first
    best = get(refresh_best!(w; force = true), n, typemax(Int))
    npts < best || return
    t = @elapsed θr, ok, cond = refine(w.d, s, n, θ)
    wmin, mmin = margin(w.d, s, θr)
    if !(ok && wmin > 0 && mmin > 0)
        say(w, @sprintf("degree %d: %d points found (%s) but refinement failed", n, npts, how))
        return
    end
    sep = min_separation(s, θr)
    if sep <= MIN_SEPARATION                  # two orbits about to merge: a move for prune or finish, not a rule
        say(w, @sprintf("degree %d: %d points found (%s) but nodes %.1e apart", n, npts, how, sep))
        return
    end
    best = get(refresh_best!(w; force = true), n, typemax(Int))     # someone may have beaten it meanwhile
    npts < best || return
    note = best == typemax(Int) ? "found by scripts/search_campaign (pooled)" :
           "found by scripts/search_campaign (pooled); $(best - npts) point(s) fewer than the previous entry"
    e = entry(w.d, n, s, θr, cond, note)
    e["found"] = stamp(); e["worker"] = w.id; e["how"] = how; e["excess"] = excess
    w.counter += 1
    atomic_write(io -> TOML.print(io, e; sorted = true),
                 joinpath(candir(w), @sprintf("deg%03d_p%06d_w%d_%d.toml", n, npts, w.id, w.counter)))
    w.best[n] = npts
    w.excess[n] = excess
    write_best_table(w)
    say(w, best == typemax(Int) ? @sprintf("NEW degree %d: %d points (excess %d, %s, refined in %.0f s)", n, npts, excess, how, t) :
                                  @sprintf("IMPROVED degree %d: %d → %d points (excess %d, %s)", n, best, npts, excess, how))
end

# --- a positive-weight solve: weights as wᵢ = w̄ᵢ exp(xᵢ) ------------------------------------------

"""
Levenberg–Marquardt with every weight kept positive through `w = w̄ exp(x)`, steps in the
log-weights capped at 2, and a stagnation stop: after `window` accepted steps without a 2%
improvement above the noise floor it gives up. Returns (θ, residual).
"""
function positive_solve(d, s, θ0, n; maxiter = 300, tol = 1e-13, window = 40)
    sys = system(d, s, n)
    widx = weight_indices(d, s)
    wbar = max.(abs.(θ0[widx]), 1e-14)
    x = copy(θ0)
    x[widx] .= log.(max.(θ0[widx], 1e-3 .* wbar) ./ wbar)
    θof(x) = (θ = copy(x); θ[widx] .= wbar .* exp.(x[widx]); θ)
    function F(x)
        θ = θof(x)
        r, J = sys(θ)
        J[:, widx] .*= transpose(θ[widx])
        return r, J
    end
    r, J = F(x)
    nr = norm(r)
    μ = 1e-3 * maximum(sum(abs2, J; dims = 1))
    hist = [nr]
    for _ in 1:maxiter
        nr <= tol && break
        A = J' * J
        F_ = lu(A + μ * Diagonal(diag(A) .+ 1e-12); check = false)
        if !issuccess(F_)
            μ *= 8; μ > 1e20 && break; continue
        end
        Δ = F_ \ (J' * r)
        Δ .*= min(1.0, 2.0 / max(maximum(abs, view(Δ, widx)), 1e-12))
        xn = x - Δ
        if !all(isfinite, xn) || !inside(d, s, θof(xn))
            μ *= 8; μ > 1e20 && break; continue
        end
        rn, Jn = F(xn)
        nrn = norm(rn)
        if nrn < nr
            x, r, J, nr = xn, rn, Jn, nrn
            μ = max(μ / 3, 1e-15)
            push!(hist, nr)
            k = length(hist)
            (k > window && nr > 1e-8 && nr > 0.98 * hist[k - window]) && break
        else
            μ *= 4; μ > 1e20 && break
        end
    end
    return θof(x), nr
end

# --- the three roles -----------------------------------------------------------------------------

"Grow: from the best rule one degree step down, add orbits to `E` unknowns above square, and fit."
function grow!(w, n)
    d = w.d
    base = get(best_entries(w), n - degree_step(d), nothing)
    base === nothing && return false
    s0, θ0 = structure(d, base["structure"]), Float64.(base["seed"])
    m = n_equations(d, n)
    E = grow_excess(m, w.rng)
    added = Any[]
    u = CR.nunknowns(s0)
    if rand(w.rng) < 0.3                                   # recipe: one of each smaller type as well
        for t in small_orbits(d)
            push!(added, t); u += CR.nunknowns(t)
        end
    end
    g = general_orbit(d)
    while u < m + E
        push!(added, g); u += CR.nunknowns(g)
    end
    for t in 1:8
        s, θ = with_orbits(d, s0, θ0, [random_orbit(d, o, w.rng) for o in added])
        r = fit(d, s, θ, n)
        r === nothing && continue
        w.stats[:grow] += 1
        say(w, @sprintf("grow degree %d: %d points, excess %d (try %d)", n, CR.npoints(r[1]), CR.nunknowns(r[1]) - m, t))
        record!(w, n, r[1], r[2], "grow")
        return true
    end
    say(w, "grow degree $n: no fit from $(length(added)) added orbits, target excess $E")
    return false
end

"""
Chain: one of the package's own grow → eliminate chains from the best rule a degree step down.
Where a chain costs seconds (low degrees, boxes) these are hard to beat; the result is kept like
any other find.
"""
function chain!(w, n)
    base = get(best_entries(w), n - degree_step(w.d), nothing)
    base === nothing && return false
    t = @elapsed r = chain(w.d, structure(w.d, base["structure"]), Float64.(base["seed"]), n, rand(w.rng, UInt32))
    w.stats[:chain] = get(w.stats, :chain, 0) + 1
    r === nothing && (say(w, @sprintf("chain degree %d: nothing (%.0f s)", n, t)); return false)
    say(w, @sprintf("chain degree %d: %d points, excess %d (%.0f s)", n, CR.npoints(r[1]),
                    CR.nunknowns(r[1]) - n_equations(w.d, n), t))
    record!(w, n, r[1], r[2], "chain")
    return true
end

"""
Prune: from a parent far above square, remove points step by step as the package's
`eliminate` does — moves removing most points first, in random order among those — refitting
after each and saving every step, as far as any move succeeds. The rule where it stops is
marked exhausted: further moves from it are left to the finish step, which perturbs them.
"""
function prune!(w, n, parent)
    d = w.d
    _, s, θ, _ = load_rule(w, parent.path)
    m = n_equations(d, n)
    steps = 0
    while CR.nunknowns(s) > m
        mvs = filter(mv -> CR.nunknowns(mv[1]) >= m, moves(d, s, θ))
        isempty(mvs) && break
        sort!(mvs; by = mv -> (-mv[3], rand(w.rng)))      # a different path for every worker
        ok = false
        for mv in mvs[1:min(PRUNE_TRIES, length(mvs))]
            r = fit(d, mv[1], mv[2], n)
            r === nothing && continue
            s, θ = r
            record!(w, n, s, θ, "prune")
            steps += 1; ok = true
            break
        end
        ok || (push!(w.exhausted, pool_path(w, n, s, θ)); break)
    end
    w.stats[:prune] += steps
    say(w, @sprintf("prune degree %d: %d steps from %d to %d points, excess %d", n, steps, parent.npoints,
                    CR.npoints(s), CR.nunknowns(s) - m))
end

"Finish: one cheap attempt at removing points from a near-square parent, chosen family by family."
function finish!(w, n, near)
    d = w.d
    fam(p) = get!(() -> family(d, load_rule(w, p.path)[2]), w.families, p.path)
    groups = Dict{Any,Vector{Any}}()
    for p in near
        push!(get!(groups, fam(p), Any[]), p)
    end
    keys_ = collect(keys(groups))
    counts = [get(w.family_tries, k, 0) for k in keys_]
    key = rand(w.rng, keys_[counts .== minimum(counts)])
    w.family_tries[key] = get(w.family_tries, key, 0) + 1
    members = sort(groups[key]; by = p -> (p.npoints, get(w.tries, p.path, 0)))
    parent = members[rand(w.rng, 1:min(5, length(members)))]
    k = (w.tries[parent.path] = get(w.tries, parent.path, 0) + 1)
    _, s, θ, _ = load_rule(w, parent.path)
    m = n_equations(d, n)
    mvs = filter(mv -> CR.nunknowns(mv[1]) >= m, moves(d, s, θ))
    isempty(mvs) && return
    most = maximum(mv -> mv[3], mvs)
    pool = rand(w.rng) < 0.7 ? filter(mv -> mv[3] == most, mvs) : mvs
    mv = rand(w.rng, pool)
    θ0 = copy(mv[2])
    σ = SIGMAS[mod1(k, length(SIGMAS))]
    if σ > 0
        widx = Set(weight_indices(d, mv[1]))
        for i in eachindex(θ0)
            i in widx || (θ0[i] += σ * randn(w.rng))
        end
    end
    w.stats[:finish] += 1
    r = if iseven(k)
        θp, nr = positive_solve(d, mv[1], θ0, n)
        nr < 1e-8 ? fit(d, mv[1], θp, n) : nothing
    else
        fit(d, mv[1], θ0, n)
    end
    if r === nothing
        w.stats[:finish_fail] += 1
    else
        w.stats[:finish_ok] += 1
        record!(w, n, r[1], r[2], iseven(k) ? "finish, positive" : "finish")
    end
end

# --- the worker loop -------------------------------------------------------------------------------

"""
The degree to work at. In extend mode, the lowest degree above the shipped table whose best rule
is still far above square (a grown rule not yet pruned), until this worker has spent `PATIENCE`
hours there without improving it; otherwise one degree step above the best table. Moving up from
an unsettled rule would grow every later degree from a poor base.
"""
function target_degree(w, mode, lo, hi, maxdeg)
    best = refresh_best!(w)
    top = maximum(keys(best))
    step = degree_step(w.d)
    if mode == "extend"
        for n in (w.shipped_top + step):step:top
            haskey(best, n) || continue
            get(w.excess, n, 0) > FINISH_EXCESS && get(w.worked, (n, best[n]), 0.0) < 3600PATIENCE && return n
        end
        top + step <= maxdeg && return top + step
    end
    lo = lo === nothing ? max(minimum(keys(best)), top - 3step) : lo
    hi = hi === nothing ? top : hi
    return rand(w.rng, lo:step:hi)
end

function run_worker(domain, root, id; mode = "extend", lo = nothing, hi = nothing, maxdeg = 1000, hours = Inf,
                    seed = 0x5ea7c4, blind = typemax(Int))
    d = DOMAINS[domain]
    dir = joinpath(root, domain)
    mkpath(joinpath(dir, "logs"))
    log = open(joinpath(dir, "logs", "worker$(id).log"), "a")
    w = Worker(d, dir, id, MersenneTwister(seed + 7919id), log, 0, Dict{String,Int}(), Dict{Any,Int}(),
               Dict{String,Any}(), Dict{Int,Int}(), Dict{Int,Int}(), Dict{Tuple{Int,Int},Float64}(), 0,
               0.0, Dict(:grow => 0, :prune => 0, :finish => 0, :finish_ok => 0, :finish_fail => 0),
               blind, Set{String}())
    w.shipped_top = maximum(e["degree"] for e in shipped(w))
    say(w, "worker $id started: $domain, mode $mode, pid $(getpid())")
    started = time()
    iter = 0
    chain_time, total_time = 0.0, 0.0
    while !isfile(joinpath(root, "STOP")) && !isfile(joinpath(dir, "STOP"))
        (time() - started) / 3600 < hours || break
        iter += 1
        n = target_degree(w, mode, lo, hi, maxdeg)
        b0 = get(w.best, n, 0)
        pl = pool_list(w, n)
        high = filter(p -> p.excess > FINISH_EXCESS && !(p.path in w.exhausted), pl)
        near = filter(p -> 0 < p.excess <= FINISH_EXCESS || p.path in w.exhausted, pl)
        roll = rand(w.rng)
        t0 = time()
        try
            if chain_time < CHAIN_SHARE * total_time       # chains get a share of the time, not of the tasks
                chain!(w, n)
                chain_time += time() - t0
            elseif isempty(pl) || (length(pl) < 3 && roll < 0.5)
                grow!(w, n)
            elseif isempty(near) || (!isempty(high) && roll < 0.2)
                isempty(high) ? grow!(w, n) : prune!(w, n, rand(w.rng, sort(high; by = p -> p.npoints)[1:min(3, length(high))]))
            elseif roll < 0.28
                grow!(w, n)
            else
                finish!(w, n, near)
            end
        catch e
            e isa InterruptException && rethrow()
            say(w, "error at degree $n: " * first(sprint(showerror, e), 300))
        end
        total_time += time() - t0
        w.worked[(n, b0)] = get(w.worked, (n, b0), 0.0) + time() - t0
        if iter % 50 == 0
            say(w, "progress: " * join(("$k=$v" for (k, v) in sort(collect(w.stats))), ", ") *
                   @sprintf(", %.2f h", (time() - started) / 3600))
        end
    end
    say(w, "worker $id stopped")
    close(log)
end
