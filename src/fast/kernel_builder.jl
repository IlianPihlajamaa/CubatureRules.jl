# Fast rules for the boundary-element kernels on flat triangles (notes/fast-kernel-rules.md):
# the constructions of DuffyGauss, DuffySinh and DuffyFinitePart in Float64, with every size
# fixed in advance instead of refined until two resolutions agree, so that a rule takes
# microseconds rather than tenths of a second. What a certified rule shows for itself, these
# show through an offline validation against the certified rules over the geometries they
# accept (scripts/validate_fast_kernel_rules.jl, src/data/fast_kernel_validation.toml); a
# geometry outside that range is handed to the certified construction.

"""
    KernelRuleBuilder(kernel, degree; dim = 3)

Builds rules for the boundary-element kernels on flat triangles fast enough to build one per
pair of collocation point and panel: `kernel` is [`InverseDistance`](@ref),
[`InverseDistanceCubed`](@ref) or [`InverseDistanceGradient`](@ref) (the type), `degree` the
polynomial degree, `dim` 3 for triangles in space ([`SurfaceTriangle`](@ref)) or 2 for
triangles in the plane. Only `Float64` rules, degrees 0 to 19.

```julia
b = KernelRuleBuilder(InverseDistance, 9)
T = SurfaceTriangle((0, 0, 0), (1, 0, 0), (0, 1, 0))
x, w = b(T, InverseDistance((0.3, 0.3, 1e-4)))      # views into b's buffers
integrate(y -> 1 + y[1], b, T, InverseDistance((0.3, 0.3, 0.0)))
```

The rules are the ones the certified families build — [`DuffyGauss`](@ref) and
[`DuffyFinitePart`](@ref) for `x₀` on the triangle, [`DuffySinh`](@ref) off it — with the same
point counts and degree claims, but built in Float64 with their discretisations sized in
advance, in microseconds. They are not verified one by one: the builder was validated
against the certified rules over a sweep of triangles and points (see `validation(b)`), and
a geometry outside the validated range, such as a sub-triangle thinner than the sweep
covered, is built by the certified construction instead (`b.fallbacks` counts them).
`rule(b, T, kernel)` returns an ordinary [`QuadratureRule`](@ref) that `check` can verify.

A builder reuses its buffers: each call overwrites the nodes and weights of the last, and
one builder must not be shared between threads.
"""
mutable struct KernelRuleBuilder{K,D}
    degree::Int
    claimed::Int
    m::Int                      # points per ray (DuffyGauss, DuffySinh) or per ray direction (finite part)
    work::QuadraticWeightWork
    nodes::Vector{SVector{D,Float64}}
    weights::Vector{Float64}
    sx::Vector{Float64}         # Gauss–Legendre on [0, 1], for DuffyGauss's radial direction
    sw::Vector{Float64}
    fs::Vector{Float64}         # the finite part's radial nodes, and its weights ω and χ (see DuffyFinitePart)
    fω::Vector{Float64}
    fχ::Vector{Float64}
    fallbacks::Int
    place::Symbol
end

const FAST_KERNELS = (InverseDistance, InverseDistanceCubed, InverseDistanceGradient)

function KernelRuleBuilder(kernel::Type, degree::Integer; dim::Integer = 3)
    kernel in FAST_KERNELS || throw(ArgumentError("KernelRuleBuilder takes InverseDistance, InverseDistanceCubed or " *
                                                  "InverseDistanceGradient, got $kernel"))
    0 <= degree <= 19 || throw(ArgumentError("KernelRuleBuilder covers degrees 0 to 19, got $degree"))
    dim in (2, 3) || throw(ArgumentError("dim must be 2 or 3, got $dim"))
    fast_gl_init()
    d = Int(degree)
    gradient = kernel === InverseDistanceGradient
    msinh = gradient ? duffy_m(d + 1) : duffy_m(d)
    m = kernel === InverseDistance ? duffy_m(d) : msinh
    claimed = kernel === InverseDistance ? 2m - 1 : d
    x, w = setprecision(() -> gauss_jacobi_work(duffy_m(d), 0, 0, 128)[1:2], BigFloat, 128)
    sx, sw = Float64.((1 .+ x) ./ 2), Float64.(w ./ 2)
    fs, fω, fχ = kernel === InverseDistance ? (Float64[], Float64[], Float64[]) : finitepart_radial_weights(kernel, d)
    mmax = max(msinh, m, cld(d + 2, 2), 2) + 1
    return KernelRuleBuilder{kernel,dim}(d, claimed, m, QuadraticWeightWork(mmax), SVector{dim,Float64}[], Float64[],
                                         sx, sw, fs, fω, fχ, 0, :none)
end

Base.show(io::IO, b::KernelRuleBuilder{K,D}) where {K,D} =
    print(io, "KernelRuleBuilder(", K, ", ", b.degree, "; dim = ", D, ")")

# The finite part's radial rule on degree + 1 Gauss–Legendre nodes: ω for the finite part in s
# of the polynomial part, χ for φ(x₀) or ∇φ(x₀)·v from the values along a ray (as in
# DuffyFinitePart's build, in BigFloat and rounded).
function finitepart_radial_weights(kernel, d)
    hyper = kernel === InverseDistanceCubed
    ms = d + 1
    bits = 53 + 8d + 48
    return setprecision(BigFloat, bits) do
        sx, _ = gauss_jacobi_work(ms, 0, 0, bits)[1:2]
        s = (1 .+ sx) ./ 2
        Vs = [s[i]^k for k in 0:(ms - 1), i in 1:ms]
        weights_for(μ) = Vs \ [μ(k) for k in 0:(ms - 1)]
        ω = hyper ? weights_for(k -> k == 0 ? -one(BigFloat) : k == 1 ? zero(BigFloat) : 1 / BigFloat(k - 1)) :
                    weights_for(k -> k == 0 ? zero(BigFloat) : 1 / BigFloat(k))
        χ = weights_for(k -> k == (hyper ? 1 : 0) ? one(BigFloat) : zero(BigFloat))
        Float64.(s), Float64.(ω), Float64.(χ)
    end
end

# ---- geometry ----------------------------------------------------------------------------------

_fcross(a::SVector{2}, b::SVector{2}) = abs(a[1] * b[2] - a[2] * b[1])
_fcross(a::SVector{3}, b::SVector{3}) = norm(cross(a, b))
# twice the signed area of (a, b) along the normal n (in the plane: n = nothing)
_fnormal(a::SVector{2}, b::SVector{2}) = nothing
_fnormal(a::SVector{3}, b::SVector{3}) = cross(a, b)
_fsarea(a::SVector{2}, b::SVector{2}, ::Nothing) = a[1] * b[2] - a[2] * b[1]
_fsarea(a::SVector{3}, b::SVector{3}, n) = dot(n, cross(a, b))
# the point projected onto the plane, and whether it was within `tol` of it
_fproject(X::SVector{2}, V1, ::Nothing, tol) = (X, true)
function _fproject(X::SVector{3}, V1, n, tol)
    nn = dot(n, n)
    σ = dot(n, X - V1)
    return X - (σ / nn) * n, σ^2 <= tol^2 * nn
end

const FAST_EDGES = ((2, 3), (3, 1), (1, 2))

# a piece's edge vectors from the apex x, rounded once (zero for an unused piece)
function _fpiece(V::SVector{3,SVector{D,T}}, x::SVector{D,T}, e::Int) where {D,T}
    e == 0 && return (zero(SVector{D,Float64}), zero(SVector{D,Float64}))
    j, k = FAST_EDGES[e]
    return SVector{D,Float64}(V[j] - x), SVector{D,Float64}(V[k] - x)
end

"""
    fast_kernel_geometry(V, X) -> (; place, c, w, npatches, edges, pieces, moved)

[`kernel_geometry`](@ref) for the vertices `V` and the point `X`, in Float64x2: `place` is
`:vertex`, `:edge` or `:interior` for `X` on the triangle (within 64 ulps of its plane and
boundary, and then moved there by `moved`) and `:near` off it; `c` is the point the triangle is
cut at, rounded to Float64, and `w ≈ c − X` its offset from `X`, accurate to a few ulps of
itself however close `X` is (`n·(X − V₁)` in Float64 would be wrong by `eps |X|`, which is
everything for `1/r³` at a height of `1e-12`); the pieces are the edges
`FAST_EDGES[edges[i]]`, `i ≤ npatches`, not through `c`, and `pieces[i]` their edge vectors
from the unrounded apex, rounded once.
"""
function fast_kernel_geometry(V64::SVector{3,SVector{D,Float64}}, X64::SVector{D,Float64}) where {D}
    tol = 64 * max(maximum(x -> eps(abs(x)), X64), maximum(v -> maximum(x -> eps(abs(x)), v), V64))
    tol2 = Float64x2(tol)^2
    V = map(v -> SVector{D,Float64x2}(v), V64)
    X = SVector{D,Float64x2}(X64)
    n = _fnormal(V[2] - V[1], V[3] - V[1])
    total = _fsarea(V[2] - V[1], V[3] - V[1], n)
    Xp, inplane = _fproject(X, V[1], n, tol)
    zero2, one2 = zero(Float64x2), one(Float64x2)
    inside = true
    best, c, ce, ct = Float64x2(Inf), Xp, 0, zero2
    for (i, (j, k)) in enumerate(FAST_EDGES)
        _fsarea(V[j] - Xp, V[k] - Xp, n) * total < 0 && (inside = false)
        d = V[k] - V[j]
        t = clamp(dot(Xp - V[j], d) / dot(d, d), zero2, one2)
        p = V[j] + t * d
        q = dot(Xp - p, Xp - p)
        if q < best
            best, c, ce, ct = q, p, i, t
        end
    end
    # where the cut point lies: at a vertex, on an edge, or inside
    vertex, edge = 0, 0
    x = Xp
    if best <= tol2 || !inside
        x = c
        edge = ce
        ct == zero2 && (vertex = FAST_EDGES[ce][1])
        ct == one2 && (vertex = FAST_EDGES[ce][2])
        if best <= tol2
            for i in 1:3
                if dot(x - V[i], x - V[i]) <= tol2
                    x, vertex = V[i], i
                end
            end
        end
    end
    edges = vertex > 0 ? (vertex, 0, 0) : edge > 0 ? (mod1(edge + 1, 3), mod1(edge + 2, 3), 0) : (1, 2, 3)
    npatches = vertex > 0 ? 1 : edge > 0 ? 2 : 3
    cf = SVector{D,Float64}(x)
    w = SVector{D,Float64}(x - X)
    # each piece's edge vectors from the unrounded apex: rounded only once, relative to
    # themselves, so the kernel sees x₀ moved by eps |a| rather than eps |X|
    pieces = (_fpiece(V, x, edges[1]), _fpiece(V, x, edges[2]), _fpiece(V, x, edges[3]))
    if inplane && (inside || best <= tol2)
        place = npatches == 1 ? :vertex : npatches == 2 ? :edge : :interior
        return (; place, c = cf, w, npatches, edges, pieces, moved = Float64(norm(x - X)))
    end
    return (; place = :near, c = cf, w, npatches, edges, pieces, moved = 0.0)
end

# ---- the three constructions -------------------------------------------------------------------

# the validated range: below these, a piece goes to the certified construction
const FAST_MIN_ETA = 1e-8               # DuffyGauss: the piece's η = 2 area / |d|²
const FAST_MIN_ETA_FINITEPART = 1e-3    # DuffyFinitePart
const FAST_MIN_KAPPA_SINH = 1e-6        # DuffySinh: x₀'s distance from the far edge's line, over |d|

"""
    sinh_angular_lines(K, L, strength = 1) -> n

Gauss–Legendre lines per piece for DuffySinh's angular direction, exact to degree `K`, `L`
half the length of the interval in `v`: the integrand is analytic in the strip
`|Im v| < π/2`, so Gauss–Legendre converges like `ρ^{-2n}` with `ρ` the largest Bernstein
ellipse in it; calibrated with a margin against the certified rules. The degree term has a
floor: at degree 0 the gradient kernel, whose angular integrand behaves like `q^{-3/2}` with a
logarithm and cancels as a principal value for `e` in the plane, needed four lines more than
the formula without it gave.

The singularity on the strip's edge is the far edge's point nearest `x₀`. For `1/r` it is
logarithmic; for `1/r³` and its gradient it is `strength = η/κ²` times larger than the
integral (`η` the piece's in-plane distance to the far edge, `κ` the distance in space, both
over the edge's length), and the error grows with it. This term matters for `x₀` just above
the inside of a thin triangle, close to two long edges: there η ≈ κ ≈ 0.006, and without it
the gradient at degree 0 was off by 6e-12 where the rounded certified rule is right to 3e-16.
"""
function sinh_angular_lines(K, L, strength = 1.0)
    b = π / (2L)
    ρ = b + sqrt(1 + b^2)
    return fast_gl_size(ceil(Int, 4 + max(K, 6) / 1.5 + 0.85 * log(max(strength, 1.0) / 4e-15) / (2log(ρ))))
end

function fast_duffygauss!(b::KernelRuleBuilder{K,D}, g) where {K,D}
    W, m = b.work, b.m
    x0 = g.c
    for i in 1:(g.npatches)
        a, bb = g.pieces[i]
        d = bb - a
        dd = dot(d, d)
        J = _fcross(a, bb)
        η = J / dd
        η >= FAST_MIN_ETA || return false
        quadratic_weight_rule!(W, -dot(a, d) / dd, η, m; j = 0, p = 1)
        sd = sqrt(dd)
        for jt in 1:m, is in 1:m
            push!(b.nodes, x0 + b.sx[is] * (a + W.x[jt] * d))
            push!(b.weights, b.sw[is] * W.w[jt] / sd * J)
        end
    end
    return true
end

function fast_duffysinh!(b::KernelRuleBuilder{K,D}, g, kernel) where {K,D}
    W, m = b.work, b.m
    power = K === InverseDistance ? 1 : 3
    gradient = K === InverseDistanceGradient
    edir = gradient ? SVector{D,Float64}(kernel.e) : zero(SVector{D,Float64})
    c = g.c
    w = g.w                                   # c − x₀, accurate however close x₀ is
    we = dot(w, edir)
    for i in 1:(g.npatches)
        a, bb = g.pieces[i]
        d = bb - a
        dd = dot(d, d)
        J = _fcross(a, bb)
        # t = t₁ + κ sinh v about the complex points where the edge point meets x₀
        aw = a + w
        ctr = -dot(aw, d) / dd
        κ = _fcross(d, aw) / dd
        κ >= FAST_MIN_KAPPA_SINH || return false
        v0, v1 = asinh(-ctr / κ), asinh((1 - ctr) / κ)
        vm, hv = (v0 + v1) / 2, (v1 - v0) / 2
        nt = sinh_angular_lines(b.claimed, hv, power == 3 ? J / dd / κ^2 : 1.0)
        vx, _, vw = fast_gl(nt)
        for (vk, wk) in zip(vx, vw)
            sv = vm + hv * vk
            t = ctr + κ * sinh(sv)
            dtdv = hv * κ * cosh(sv)
            e = a + t * d
            qq = dot(e, e)
            s0 = -dot(e, w) / qq
            κr = _fcross(e, w) / qq
            κr > 0 || s0 < 0 || return false
            quadratic_weight_rule!(W, s0, κr, m; j = 1, p = power)
            scale = wk * dtdv * J / (power == 1 ? sqrt(qq) : qq * sqrt(qq))
            ee = dot(e, edir)
            for is in 1:m
                si = W.x[is]
                push!(b.nodes, c + si * e)
                push!(b.weights, scale * W.w[is] * (gradient ? we + si * ee : 1.0))
            end
        end
    end
    return true
end

# ∫₀¹ g(t) q(t)^{-3/2} dt on the grid of q^{-3/2}'s Gauss rule, for the finite part's
# logarithmic terms; g is analytic in the same strip as DuffySinh's angular integrand
function fast_inverse_cube_integral(W, ts, η, dd, g, m)
    L = quadratic_weight_halflength(ts, η)
    bb = π / (2L)
    ρ = bb + sqrt(1 + bb^2)
    M = max(quadratic_weight_gridsize(m, L, 0, 3), fast_gl_size(ceil(Int, 12 + 0.85 * log(1 / 4e-15) / (2log(ρ)))))
    quadratic_weight_grid!(W, M, ts, η, 0, 3)
    total = g(W.s[1]) * W.λ[1]
    for i in 2:M
        total += g(W.s[i]) * W.λ[i]
    end
    return total / (dd * sqrt(dd))
end

function fast_finitepart!(b::KernelRuleBuilder{K,D}, g, kernel) where {K,D}
    W = b.work
    hyper = K === InverseDistanceCubed
    d = b.degree
    mt = hyper ? max(2, cld(d + 1, 2)) : max(1, cld(d + 2, 2))
    e = hyper ? zero(SVector{D,Float64}) : SVector{D,Float64}(kernel.e)
    x0 = g.c
    for i in 1:(g.npatches)
        A, B = g.pieces[i]
        Dv = B - A
        dd = dot(Dv, Dv)
        J = _fcross(A, B)
        η = J / dd
        η >= FAST_MIN_ETA_FINITEPART || return false
        ts = -dot(A, Dv) / dd
        Vt(t) = A + t * Dv
        # the logarithmic term, before the angular rule overwrites the grid
        Lv = hyper ? fast_inverse_cube_integral(W, ts, η, dd, t -> Vt(t) * log(norm(Vt(t))), mt) :
                     fast_inverse_cube_integral(W, ts, η, dd, t -> dot(Vt(t), e) * log(norm(Vt(t))), mt)
        quadratic_weight_rule!(W, ts, η, mt; j = 0, p = 3)
        λs = dd * sqrt(dd)
        # Σ_j κ_j v(t_j) = L, in the plane of A and B, minimum norm: κ = Vᵀ (V Vᵀ)⁻¹ L
        u1 = A / norm(A)
        u2 = B - dot(B, u1) * u1
        u2 /= norm(u2)
        y1, y2 = 0.0, 0.0
        if hyper
            G11, G12, G22 = 0.0, 0.0, 0.0
            for jt in 1:mt
                v = Vt(W.x[jt])
                p1, p2 = dot(v, u1), dot(v, u2)
                G11 += p1 * p1; G12 += p1 * p2; G22 += p2 * p2
            end
            l1, l2 = dot(Lv, u1), dot(Lv, u2)
            det = G11 * G22 - G12^2
            y1, y2 = (G22 * l1 - G12 * l2) / det, (G11 * l2 - G12 * l1) / det
        end
        sumλ = 0.0
        for jt in 1:mt
            sumλ += W.w[jt] / λs
        end
        for jt in 1:mt
            tj = W.x[jt]
            λj = W.w[jt] / λs
            v = Vt(tj)
            κj = hyper ? y1 * dot(v, u1) + y2 * dot(v, u2) : Lv * λj / sumλ
            radial = hyper ? λj : dot(v, e) * λj
            for is in eachindex(b.fs)
                push!(b.nodes, x0 + b.fs[is] * v)
                push!(b.weights, J * (radial * b.fω[is] + κj * b.fχ[is]))
            end
        end
    end
    return true
end

# ---- the call ------------------------------------------------------------------------------------

_fast_vertices(T::SurfaceTriangle) = SVector{3,SVector{3,Float64}}(map(v -> SVector{3,Float64}(v), Tuple(T.vertices)))
_fast_vertices(T::Simplex{2}) = SVector{3,SVector{2,Float64}}(map(v -> SVector{2,Float64}(v), Tuple(vertices(T))))

"""
    (b::KernelRuleBuilder)(T, kernel) -> (nodes, weights)

The rule for `WeightedDomain(T, kernel)` at `b`'s degree, into `b`'s buffers: the returned
vectors are overwritten by the next call. `T` is a [`SurfaceTriangle`](@ref) (`dim = 3`) or a
triangle in the plane (`dim = 2`); `kernel` an instance of `b`'s kernel type.
"""
function (b::KernelRuleBuilder{K,D})(T::Union{SurfaceTriangle,Simplex{2}}, kernel::PointKernel{D}) where {K,D}
    kernel isa K || throw(ArgumentError("this builder is for $K kernels, got $(typeof(kernel))"))
    V = _fast_vertices(T)
    X = SVector{D,Float64}(kernel.x0)
    g = fast_kernel_geometry(V, X)
    empty!(b.nodes)
    empty!(b.weights)
    b.place = g.place
    ok = if g.place === :near
        fast_duffysinh!(b, g, kernel)
    elseif K === InverseDistance
        fast_duffygauss!(b, g)
    else
        fast_finitepart!(b, g, kernel)
    end
    ok && return b.nodes, b.weights
    # outside the validated range: the certified construction
    b.fallbacks += 1
    r = rule(WeightedDomain(T, kernel); degree = b.degree)
    empty!(b.nodes)
    empty!(b.weights)
    append!(b.nodes, nodes(r))
    append!(b.weights, weights(r))
    b.place = :fallback
    return b.nodes, b.weights
end

"""
    integrate(f, b::KernelRuleBuilder, T, kernel)

`∫_T f(y) K(y) dy` with the rule `b` builds for `T` and `kernel`.
"""
function integrate(f, b::KernelRuleBuilder, T::Union{SurfaceTriangle,Simplex{2}}, kernel::PointKernel)
    x, w = b(T, kernel)
    total = w[1] * f(x[1])
    for i in 2:length(x)
        total += w[i] * f(x[i])
    end
    return total
end

"""
    rule(b::KernelRuleBuilder, T, kernel) -> QuadratureRule

The rule `b` builds for `T` and `kernel`, as an ordinary [`QuadratureRule`](@ref) with its
provenance and with a certificate that refers to the builder's validation instead of a
residual of its own: `check` verifies it like any other rule.
"""
function rule(b::KernelRuleBuilder{K,D}, T::Union{SurfaceTriangle,Simplex{2}}, kernel::PointKernel{D}) where {K,D}
    x, w = b(T, kernel)
    b.place === :fallback && return rule(WeightedDomain(T, kernel); degree = b.degree)
    fam = b.place === :near ? "DuffySinh" : K === InverseDistance ? "DuffyGauss" : "DuffyFinitePart"
    val = validation(b)
    worst = val === nothing ? NaN : val.worst
    prov = Provenance(family = fam * " (fast)", derivation = Derived(),
                      path = ["KernelRuleBuilder: the construction of $fam in Float64, with its discretisations " *
                              "sized in advance (notes/fast-kernel-rules.md); x₀ " *
                              (b.place === :near ? "off the triangle" : "on the triangle (" * string(b.place) * ")"),
                              val === nothing ? "not validated at this degree" :
                              "validated against the certified rules: worst relative error " *
                              @sprintf("%.1e", val.worst) * " of the Dubiner moments over $(val.geometries) " *
                              "geometries ($(val.source))"],
                      seed_source = "none (derived)", citations = [DUFFY_1982, JOHNSTON_ELLIOTT_2005])
    cert = Certificate(equations = "none for this rule: the fast construction, validated offline over its range",
                       residual = BigFloat(isnan(worst) ? Inf : worst; precision = 64), residual_bits = 53,
                       digits = 15, guard_digits = 0)
    return QuadratureRule(copy(x), copy(w), WeightedDomain(T, kernel), PolynomialDegree(b.claimed), prov, cert)
end

# ---- the validation record ---------------------------------------------------------------------

const FAST_VALIDATION_FILE = joinpath(@__DIR__, "..", "data", "fast_kernel_validation.toml")
include_dependency(FAST_VALIDATION_FILE)
const FAST_VALIDATION = isfile(FAST_VALIDATION_FILE) ? TOML.parsefile(FAST_VALIDATION_FILE) : Dict{String,Any}()

"""
    validation(b::KernelRuleBuilder) -> (; worst, geometries, source) or nothing

What the offline validation found for `b`'s kernel and degree: the worst relative error of the
Dubiner moments against the certified rules, the number of geometries, and the script and
date of the run.
"""
function validation(b::KernelRuleBuilder{K,D}) where {K,D}
    entry = get(get(FAST_VALIDATION, string(nameof(K)), Dict()), string(b.degree), nothing)
    entry === nothing && return nothing
    return (; worst = Float64(entry["worst"]), geometries = Int(entry["geometries"]),
            source = string(get(FAST_VALIDATION, "source", "")))
end
