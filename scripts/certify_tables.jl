# Check the shipped Float64 seed tables so that Float64 requests can be served from them
# without refinement (see `build_symmetric` and the Lebedev `build`).
#
#     julia --project -t auto scripts/certify_tables.jl [triangle|tetrahedron|lebedev|square|cube ...]
#     julia --project -t auto scripts/certify_tables.jl triangle:51,52,53 tetrahedron:17,21,22
#
# A table name alone checks every entry; `name:degrees` only those degrees, which is what
# newly merged entries need.
#
# For every entry: refine the stored seed to 160 bits, round the result to Float64, and
# compare it with the stored seed. A stored seed that is not the correctly rounded value is
# replaced by it. Then record the defining-equation residual of the stored seed at 256 bits
# (`residual`, `residual_bits`), which the certificate of a shipped rule quotes. An entry
# without a recorded residual is refined on every request, as before.

using CubatureRules, Printf, TOML
const CR = CubatureRules

const DATA = joinpath(@__DIR__, "..", "src", "data")
# file, kind of moment system, and its N (simplex), D (box) or nothing (sphere)
const TABLES = Dict("triangle" => ("triangle_s3_seeds.toml", :simplex, 3),
                    "tetrahedron" => ("tetrahedron_s4_seeds.toml", :simplex, 4),
                    "lebedev" => ("lebedev_seeds.toml", :octahedral, nothing),
                    "square" => ("square_d4_seeds.toml", :box, 2),
                    "cube" => ("cube_oh_seeds.toml", :box, 3))
const REFINE_BITS = 160
const RESIDUAL_BITS = 256
const NOTE = "# `residual` is the defining-equation residual of the stored seed at `residual_bits`, " *
             "recorded by\n# scripts/certify_tables.jl, which also checks that each seed is the correctly " *
             "rounded Float64\n# value of the rule refined to $REFINE_BITS bits."

function refine_and_residual(e, kind, N)
    n, θ64 = e["degree"], Float64.(e["seed"])
    if kind === :box
        s = CR.BoxStructure([Int[x for x in m] for m in e["structure"]], N)
        θ, res, _ = CR.refine_box(s, n, θ64, REFINE_BITS)
        system = (θ) -> CR.BoxMomentSystem(s, n, BigFloat)(θ; jacobian = false)[1]
    elseif kind === :octahedral
        s = CR.OctahedralStructure(Symbol.(e["structure"]))
        θ, res, _ = CR.refine_octahedral(s, n, θ64, REFINE_BITS)
        system = (θ) -> CR.OctahedralMomentSystem(s, n, BigFloat)(θ)[1]
    else
        s = CR.SymmetricStructure([Int.(m) for m in e["structure"]], N)
        basis = CR.invariant_basis(N, n)
        θ, res, _ = CR.refine_symmetric(s, n, θ64, REFINE_BITS; basis)
        system = (θ) -> CR.SymmetricMomentSystem(s, n, BigFloat, basis)(θ; jacobian = false)[1]
    end
    res.converged || error("degree $n did not converge")
    θr = Float64.(θ)
    ulps = maximum(abs(a - b) / eps(max(abs(a), floatmin())) for (a, b) in zip(θr, θ64))
    resid = CR.with_bits(RESIDUAL_BITS) do
        maximum(abs, system(BigFloat.(θr)))
    end
    return θr, ulps, Float64(resid)
end

function certify(name, degrees = nothing)
    file, kind, N = TABLES[name]
    path = joinpath(DATA, file)
    header = String[]
    for l in eachline(path)
        startswith(l, "#") || break
        push!(header, l)
    end
    entries = TOML.parsefile(path)["rule"]
    todo = [i for (i, e) in enumerate(entries)
            if get(e, "status", "ok") == "ok" && (degrees === nothing || e["degree"] in degrees)]
    results = Vector{Any}(undef, length(entries))
    Threads.@threads :dynamic for i in todo
        t = @elapsed results[i] = refine_and_residual(entries[i], kind, N)
        e = entries[i]
        θr, ulps, resid = results[i]
        @printf("%-11s degree %3d, %4d points: seed %s (%.1f ulp), residual %.1e  [%.1f s]\n", name, e["degree"],
                e["npoints"], θr == Float64.(e["seed"]) ? "correctly rounded" : "REPLACED", ulps, resid, t)
        flush(stdout)
    end
    for i in todo
        θr, _, resid = results[i]
        entries[i]["seed"] = θr
        entries[i]["residual"] = resid
        entries[i]["residual_bits"] = RESIDUAL_BITS
    end
    occursin(NOTE, join(header, "\n")) || append!(header, ["#", split(NOTE, "\n")...])
    open(path, "w") do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for arg in (isempty(ARGS) ? ["lebedev", "tetrahedron", "triangle", "square", "cube"] : ARGS)
        name, sel = occursin(":", arg) ? split(arg, ":"; limit = 2) : (arg, nothing)
        certify(String(name), sel === nothing ? nothing : parse.(Int, split(sel, ",")))
    end
end
