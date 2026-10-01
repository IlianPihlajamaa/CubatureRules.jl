# Import high-precision O_h sphere rules computed outside this package into the shipped seed
# table src/data/lebedev_seeds.toml, as Float64 orbit parameters.
#
#     julia --project scripts/import_octahedral_rules.jl DIR [degree ...]
#
# DIR holds, for each degree d, `degree<d>_parameters.dat` (one decimal per line, at least 100
# significant digits) and `degree<d>_metadata.json` with `"structure": [a1, a2, a3, N1, N2, N3]`.
# That is the layout of the catalogue the degree-133 to 201 entries came from. The parameter
# vector is
#
#     orbit masses, summing to 1: the enabled fixed orbits a1, a2, a3, then N1 b, N2 c, N3 d
#     N1 b polar angles θ    representative (sin θ/√2, sin θ/√2, cos θ)
#     N2 c azimuths φ        representative (cos φ, sin φ, 0)
#     N3 d polar angles θ, then N3 d azimuths φ
#                            representative (sin θ cos φ, sin θ sin φ, cos θ)
#
# Each rule is converted at 768 bits to this package's θ = [w₁, p₁…, w₂, p₂…, …], with per-point
# weights summing to 4π, then rounded to Float64. Every orbit is invariant under signed
# permutations, so the representative is free: for c and d it is taken with the coordinate
# recovered through the square root the largest (at least 1/√2 and 1/√3), which keeps that
# coordinate accurate when the seed is expanded in Float64.
#
# The entries are written without a residual. scripts/certify_tables.jl then refines each seed
# to 160 bits, checks that it is the correctly rounded value, and records the residual; until
# then a Float64 request refines instead of shipping the seed as stored.

using CubatureRules, Printf, TOML
const CR = CubatureRules

const SEED_FILE = joinpath(@__DIR__, "..", "src", "data", "lebedev_seeds.toml")

function read_structure(path)
    m = match(r"\"structure\"\s*:\s*\[([^\]]*)\]", read(path, String))
    m === nothing && error("no structure in $path")
    return parse.(Int, strip.(split(m[1], ",")))
end

"This package's orbit list and θ for the rule of `degree` stored in `dir`, rounded to `T` (768-bit BigFloat for `T = BigFloat`)."
function convert_rule(dir, degree; T = Float64)
    a1, a2, a3, N1, N2, N3 = read_structure(joinpath(dir, "degree$(degree)_metadata.json"))
    setprecision(BigFloat, 768) do
        u = [parse(BigFloat, strip(l)) for l in eachline(joinpath(dir, "degree$(degree)_parameters.dat"))
             if !isempty(strip(l))]
        nfix = a1 + a2 + a3
        norb = nfix + N1 + N2 + N3
        length(u) == norb + N1 + N2 + 2N3 || error("degree $degree: $(length(u)) parameters for structure $((a1, a2, a3, N1, N2, N3))")
        masses = u[1:norb]
        θb = u[norb+1:norb+N1]
        φc = u[norb+N1+1:norb+N1+N2]
        θd = u[norb+N1+N2+1:norb+N1+N2+N3]
        φd = u[norb+N1+N2+N3+1:end]
        fourpi = 4 * BigFloat(π)
        kinds = Symbol[]
        θ = BigFloat[]
        fixed = [k for (k, on) in zip((:a1, :a2, :a3), (a1, a2, a3)) if on == 1]
        sizes = Dict(:a1 => 6, :a2 => 12, :a3 => 8, :b => 24, :c => 24, :d => 48)
        weight(i, k) = masses[i] * fourpi / sizes[k]
        for (i, k) in enumerate(fixed)
            push!(kinds, k); push!(θ, weight(i, k))
        end
        for j in 1:N1                                     # (l, l, m): l = |sin θ|/√2
            push!(kinds, :b); push!(θ, weight(nfix + j, :b), abs(sin(θb[j])) / sqrt(BigFloat(2)))
        end
        for j in 1:N2                                     # (q, √(1-q²), 0), q the smaller
            push!(kinds, :c)
            push!(θ, weight(nfix + N1 + j, :c), min(abs(cos(φc[j])), abs(sin(φc[j]))))
        end
        for j in 1:N3                                     # (r, s, √(1-r²-s²)), r ≥ s the smaller two
            x = sort(abs.([sin(θd[j]) * cos(φd[j]), sin(θd[j]) * sin(φd[j]), cos(θd[j])]))
            push!(kinds, :d); push!(θ, weight(nfix + N1 + N2 + j, :d), x[2], x[1])
        end
        abs(sum(masses) - 1) < BigFloat(10)^-60 || error("degree $degree: masses sum to $(Float64(sum(masses)))")
        return CR.OctahedralStructure(kinds), T.(θ)
    end
end

function table_entry(degree, s, θ64)
    wmin, dmin = CR.octahedral_margins(s, θ64)
    (wmin > 0 && dmin > 0) || error("degree $degree: not admissible (min weight $wmin, margin $dmin)")
    sys64 = CR.OctahedralMomentSystem(s, degree, Float64)
    r64, J64 = sys64(θ64)
    length(r64) == CR.nunknowns(s) || error("degree $degree: $(length(r64)) equations, $(CR.nunknowns(s)) unknowns")
    κ = CR.seed_cond(CR.lsq_step(J64, r64; rank_rtol = 1e-14)[2],
                     () -> CR.OctahedralMomentSystem(s, degree, BigFloat)(BigFloat.(θ64))[2])
    return Dict{String,Any}(
        "degree" => degree, "npoints" => CR.npoints(s), "structure" => String.([o.kind for o in s.orbits]),
        "seed" => θ64, "status" => "ok", "method" => "continuation", "cond" => κ,
        "min_weight" => wmin, "min_margin" => dmin,
        "note" => "square system (as many unknowns as invariant conditions); no published rule at this degree")
end

function import_rules(dir, degrees)
    header = String[]
    for l in eachline(SEED_FILE)
        startswith(l, "#") || break
        push!(header, l)
    end
    entries = TOML.parsefile(SEED_FILE)["rule"]
    for d in degrees
        t = @elapsed (s, θ64) = convert_rule(dir, d)
        e = table_entry(d, s, θ64)
        filter!(x -> x["degree"] != d, entries)
        push!(entries, e)
        @printf("degree %3d: %5d points, %3d orbits, %3d unknowns, min weight %.3e, margin %.3e, cond %.1e  [%.1f s]\n",
                d, e["npoints"], length(s.orbits), length(θ64), e["min_weight"], e["min_margin"], e["cond"], t)
        flush(stdout)
    end
    sort!(entries; by = x -> x["degree"])
    open(SEED_FILE, "w") do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", SEED_FILE)
end

"""
    next_error(dir, degree) -> Float64

How far the rule misses degree `degree + 1`, measured by `check` itself — its sharpness test
functions — on the high-precision rule in `dir`. In `Float64` these rules come too close to
the next degree for the check to resolve: the degree-`d+1` monomials are almost entirely of
lower degree on the sphere. Recorded in the table, it lets a shipped rule report that rather
than look understated (see `next_error` in `Certificate`).
"""
function next_error(dir, degree)
    s, θ = convert_rule(dir, degree; T = BigFloat)
    # 300 bits: verification holds a rule to its own precision, and the rules are accurate
    # to about 1e-100 (residual 1e-98 to 1e-114), well inside 2^-300 ≈ 5e-91
    setprecision(BigFloat, 300) do
        xs, ws = CR.expand(s, BigFloat.(θ))
        prov = CR.Provenance(family = "Lebedev", derivation = CR.Seeded(), path = String[], seed_source = dir,
                             symmetry = :Oh, license = "")
        r = CR.QuadratureRule([CR.SVector{3,BigFloat}(x) for x in xs], ws, Sphere{3}(), PolynomialDegree(degree), prov)
        v = CR.verify(r, PolynomialDegree(degree))
        (v.exact && v.sharp === true) || error("degree $degree: not exact and sharp at high precision ($(v.exact), $(v.sharp))")
        return Float64(v.sharp_residual)
    end
end

"Record `next_error` in the existing entries for `degrees`, leaving everything else as it is."
function record_next_error(dir, degrees)
    header = String[]
    for l in eachline(SEED_FILE)
        startswith(l, "#") || break
        push!(header, l)
    end
    entries = TOML.parsefile(SEED_FILE)["rule"]
    for e in entries
        e["degree"] in degrees || continue
        t = @elapsed e["next_error"] = next_error(dir, e["degree"])
        @printf("degree %3d: misses degree %d by %.2e  [%.1f s]\n", e["degree"], e["degree"] + 1, e["next_error"], t)
        flush(stdout)
    end
    open(SEED_FILE, "w") do io
        foreach(l -> println(io, l), header)
        println(io)
        TOML.print(io, Dict("rule" => entries); sorted = true)
    end
    println("wrote ", SEED_FILE)
end

if abspath(PROGRAM_FILE) == @__FILE__
    dir = ARGS[1]
    mode = "--next-error" in ARGS ? :next_error : :import
    rest = filter(!=("--next-error"), ARGS[2:end])
    degrees = !isempty(rest) ? parse.(Int, rest) :
              sort([parse(Int, m[1]) for f in readdir(dir) for m in (match(r"^degree(\d+)_parameters\.dat$", f),) if m !== nothing])
    mode === :import ? import_rules(dir, degrees) : record_next_error(dir, degrees)
end
