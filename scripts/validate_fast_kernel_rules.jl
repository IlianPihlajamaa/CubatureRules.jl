# The offline validation of KernelRuleBuilder (notes/fast-kernel-rules.md): its rules against
# the certified ones, over a sweep of triangles and points, for every kernel and degree. Writes
# src/data/fast_kernel_validation.toml, which `validation(b)` and the certificates of the fast
# rules quote.
#
#     julia --project scripts/validate_fast_kernel_rules.jl            # the full sweep (hours)
#     julia --project scripts/validate_fast_kernel_rules.jl quick      # a few geometries, two degrees
#     julia --project scripts/validate_fast_kernel_rules.jl resume     # keep the entries already recorded
#
# Triangles are drawn in normalised form — (0, 0), (1, 0) and an apex (a, b), b from 0.02 to 1.6,
# from thin to obtuse — then rotated, scaled and shifted at random in space. Points: on the
# triangle (at a vertex, on an edge, inside, some inside within 1e-7 to 1e-1 of an edge), and
# off it (above or below the triangle or beside it, at heights from 1e-12 to 3 diameters, and
# in its plane outside). For each, the fast rule and the certified rule at 20 digits integrate
# the triangle's orthonormal Dubiner polynomials to the degree; the error is the largest
# difference over the largest integral (or 1). The certified rule rounded to Float64 gives the
# floor that any Float64 rule faces, recorded beside it. Geometries the builder hands to the
# certified construction (outside its range) are counted and left out, and so are those where
# the certified construction itself fails (`unbuilt`: a piece thinner than its discretisation
# reaches, as for the finite parts with x₀ within 1e-7 of an edge).

using CubatureRules, LinearAlgebra, Printf, Random, StaticArrays, TOML
const CR = CubatureRules

const QUICK = "quick" in ARGS
const OUT = joinpath(@__DIR__, "..", "src", "data", "fast_kernel_validation.toml")
const DEGREES = QUICK ? [0, 1, 3, 9] : collect(0:19)
const N_ON, N_NEAR = QUICK ? (3, 4) : (24, 48)

function random_rotation(rng)
    q = normalize(randn(rng, 4))
    a, b, c, d = q
    return @SMatrix [a^2 + b^2 - c^2 - d^2 2(b * c - a * d) 2(b * d + a * c);
                     2(b * c + a * d) a^2 - b^2 + c^2 - d^2 2(c * d - a * b);
                     2(b * d - a * c) 2(c * d + a * b) a^2 - b^2 - c^2 + d^2]
end

function random_triangle(rng)
    a, b = -0.5 + 2rand(rng), 10^(-1.7 + 1.9rand(rng))
    P = (SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(a, b, 0.0))
    Q, scale, shift = random_rotation(rng), 10^(-1 + 2rand(rng)), SVector{3}(randn(rng, 3))
    V = map(p -> scale * (Q * p) + shift, P)
    return SurfaceTriangle(Tuple.(V)...), V
end

# a point on the triangle: a vertex, on an edge, inside, or inside near an edge
function point_on(rng, V, i)
    kind = (:vertex, :edge, :interior, :nearedge)[mod1(i, 4)]
    kind === :vertex && return V[rand(rng, 1:3)], kind
    if kind === :edge
        j = rand(rng, 1:3)
        return V[j] + (0.001 + 0.998rand(rng)) * (V[mod1(j + 1, 3)] - V[j]), kind
    end
    λ = kind === :interior ? normalize(rand(rng, 3) .+ 0.05, 1) :
        (δ = 10^(-7 + 6rand(rng)); l = [δ, rand(rng), rand(rng)]; l[2:3] .*= (1 - δ) / (l[2] + l[3]); l[randperm(rng, 3)])
    return λ[1] * V[1] + λ[2] * V[2] + λ[3] * V[3], kind
end

# a point off it: above or below an interior point at a height from 1e-12 to 3 diameters; beside
# the triangle, off an edge by 1e-6 to 0.3 diameters, at such a height; or so beside it in its plane
function point_near(rng, V, i)
    diam = maximum(norm(V[j] - V[k]) for (j, k) in ((1, 2), (2, 3), (3, 1)))
    n = normalize(cross(V[2] - V[1], V[3] - V[1]))
    h = diam * 10^(-12 + 12.5rand(rng)) * rand(rng, (-1, 1))
    kind = (:above, :above, :beside, :inplane)[mod1(i, 4)]
    if kind === :above
        λ = normalize(rand(rng, 3) .+ 0.05, 1)
        return λ[1] * V[1] + λ[2] * V[2] + λ[3] * V[3] + h * n, kind
    end
    j = rand(rng, 1:3)
    k = mod1(j + 1, 3)
    p = V[j] + rand(rng) * (V[k] - V[j])
    u = normalize(cross(V[k] - V[j], n))
    dot(u, (V[1] + V[2] + V[3]) / 3 - p) > 0 && (u = -u)       # outward
    foot = p + diam * 10^(-6 + 5.5rand(rng)) * u
    return kind === :inplane ? foot : foot + h * n, kind
end

kernel_of(K, x0, rng) = K === InverseDistance ? InverseDistance(Tuple(x0)) :
                        K === InverseDistanceCubed ? InverseDistanceCubed(Tuple(x0)) :
                        InverseDistanceGradient(Tuple(x0), Tuple(normalize(randn(rng, 3))))

# Dubiner moments, in BigFloat, of the rule (xs, ws) on T, and of the same rule rounded to Float64
function moments(V, xs, ws, K; rounded = false)
    setprecision(BigFloat, 106) do
        Vb = [SVector{3,BigFloat}(big.(v)) for v in V]
        E = hcat(Vb[2] - Vb[1], Vb[3] - Vb[1])
        A = inv(E' * E) * E'
        dub = CR.DubinerBasis{BigFloat}(K; normalize = true)
        total = zeros(BigFloat, CR.dubiner_length(K))
        for (x, w) in zip(xs, ws)
            xb = rounded ? SVector{3,BigFloat}(big.(Float64.(x))) : SVector{3,BigFloat}(big.(x))
            wb = rounded ? big(Float64(w)) : big(w)
            ξ = A * (xb - Vb[1])
            CR.dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, ξ[1], ξ[2])
            total .+= wb .* dub.φ
        end
        total
    end
end

relerr(m, ref) = Float64(maximum(abs, m - ref) / max(1, maximum(abs, ref)))

# the record, rewritten after every kernel and degree so that an interrupted run keeps what it did
function write_record(results)
    results["source"] = "scripts/validate_fast_kernel_rules.jl$(QUICK ? " quick" : ""), $(Libc.strftime("%Y-%m-%d", time())), " *
                        "$(N_ON) points on and $(N_NEAR) off random triangles per kernel and degree"
    open(OUT, "w") do io
        println(io, "# Written by scripts/validate_fast_kernel_rules.jl: the offline validation of KernelRuleBuilder")
        println(io, "# against the certified rules. Not edited by hand. `worst` is the largest relative error of the")
        println(io, "# Dubiner moments, `floor` that of the certified rule rounded to Float64, `median_us` the build time.")
        TOML.print(io, results; sorted = true)
    end
end

results = "resume" in ARGS && isfile(OUT) ? TOML.parsefile(OUT) : Dict{String,Any}()
rng = Xoshiro(20261012)                       # the triangles; the points are drawn per kernel and degree
geometries = [(random_triangle(rng)..., i) for i in 1:(N_ON + N_NEAR)]
t0 = time()
for K in (InverseDistance, InverseDistanceCubed, InverseDistanceGradient), d in DEGREES
    haskey(get(results, string(nameof(K)), Dict()), string(d)) && continue
    prng = Xoshiro(hash((string(nameof(K)), d)))
    b = KernelRuleBuilder(K, d)
    worst, floor, at, n, fallbacks, unbuilt, times = 0.0, 0.0, "", 0, 0, 0, Float64[]
    for (T, V, i) in geometries
        x0, kind = i <= N_ON ? point_on(prng, V, i) : point_near(prng, V, i)
        kern = kernel_of(K, x0, prng)
        local xs, ws, r
        try                                    # the certified construction, in a fallback, may fail too
            xs, ws = b(T, kern)
        catch err
            err isa RefinementError || rethrow()
            unbuilt += 1
            continue
        end
        if b.place === :fallback
            fallbacks += 1
            continue
        end
        xs, ws = copy(xs), copy(ws)
        push!(times, @elapsed b(T, kern))
        try
            r = setprecision(() -> rule(WeightedDomain(T, kern); degree = d, digits = 20), BigFloat, 128)
        catch err
            err isa RefinementError || rethrow()
            unbuilt += 1
            continue
        end
        ref = moments(V, nodes(r), weights(r), d)
        e = relerr(moments(V, xs, ws, d), ref)
        f = relerr(moments(V, nodes(r), weights(r), d; rounded = true), ref)
        n += 1
        floor = max(floor, f)
        if e > worst
            worst = e
            at = @sprintf("%s, x₀ %s, triangle %d", kind, b.place, i)
        end
    end
    sort!(times)
    entry = Dict("worst" => worst, "floor" => floor, "geometries" => n, "fallbacks" => fallbacks, "unbuilt" => unbuilt,
                 "worst_at" => at, "median_us" => round(1e6 * times[cld(length(times), 2)]; digits = 1))
    get!(results, string(nameof(K)), Dict{String,Any}())[string(d)] = entry
    write_record(results)
    @printf("%-24s degree %2d: worst %.1e (rounded certified %.1e) over %d geometries, %d handed on, %d unbuilt; median %.1f µs [%s] (%.0f s)\n",
            string(K), d, worst, floor, n, fallbacks, unbuilt, entry["median_us"], at, time() - t0)
    flush(stdout)
end
println("wrote ", OUT)
