# Fully symmetric rules on the square and the cube (PLAN §6 Tier 3, v0.6): orbits of the
# hyperoctahedral group B_D, the D!·2^D signed permutations of the coordinates, which is the
# symmetry group of [-1, 1]^D.
#
# An orbit is fixed by its representative's coordinates: some distinct nonzero absolute
# values, each repeated some number of times, and zeros. On the square that gives four
# types — the centre, (a, 0), (a, a) and (a, b) — and on the cube seven.
#
# The moment equations are much simpler than on the simplex. Every B_D-invariant polynomial is
# a symmetric polynomial in x₁², …, x_D², and the symmetrised products of even-degree
# orthonormal Legendre polynomials,
#
#     φ_α = c_α Σ_{σ} p_{α_σ(1)}(x₁) ⋯ p_{α_σ(D)}(x_D),     α₁ ≥ … ≥ α_D ≥ 0 all even,
#
# summed over the distinct rearrangements of α, are already orthonormal on the box. So no
# invariant basis has to be computed (compare `invariant_basis` on the simplex), only the
# constant has a nonzero integral, and the system is as well conditioned as the rule allows.
# Every orbit is centrally symmetric, so odd polynomials integrate to zero on their own and a
# rule of degree 2k is one of degree 2k + 1.

"""
    BoxOrbit(mult, D)

An orbit type of the symmetry group of `[-1, 1]^D`: the representative has `length(mult)`
distinct nonzero absolute values, the `k`th repeated `mult[k]` times, and `D - sum(mult)`
zero coordinates. `BoxOrbit([], 3)` is the centre, `BoxOrbit([2, 1], 3)` the orbit of
`(a, a, b)`.
"""
struct BoxOrbit
    mult::Vector{Int}
    D::Int
    slots::Vector{Int}      # the representative's coordinate j is parameter slots[j], 0 for zero
end
function BoxOrbit(mult::AbstractVector, D::Integer)
    m = sort(Int[x for x in mult]; rev = true)
    (all(>(0), m) && sum(m; init = 0) <= D) || throw(ArgumentError("invalid orbit $(mult) on the $(D)-cube"))
    slots = vcat(Int[], [fill(k, m[k]) for k in eachindex(m)]..., zeros(Int, D - sum(m; init = 0)))
    return BoxOrbit(m, Int(D), slots)
end
Base.:(==)(a::BoxOrbit, b::BoxOrbit) = a.D == b.D && a.mult == b.mult
Base.hash(o::BoxOrbit, h::UInt) = hash(o.mult, hash(o.D, hash(:BoxOrbit, h)))

nparams(o::BoxOrbit) = length(o.mult)
nunknowns(o::BoxOrbit) = 1 + nparams(o)
orbit_size(o::BoxOrbit) = factorial(o.D) ÷ (prod(factorial, o.mult; init = 1) * factorial(o.D - sum(o.mult; init = 0))) *
                          2^sum(o.mult; init = 0)

"Every orbit type of the symmetry group of `[-1, 1]^D`, smallest first."
function box_orbit_types(D::Integer)
    out = BoxOrbit[]
    function rec(prefix, maxpart, left)
        push!(out, BoxOrbit(prefix, D))
        for p in min(maxpart, left):-1:1
            rec(vcat(prefix, p), p, left - p)
        end
    end
    rec(Int[], D, D)
    return sort!(out; by = o -> (orbit_size(o), o.mult))
end

"""
    BoxStructure(D, orbits)

A fully symmetric rule on `[-1, 1]^D` as a list of [`BoxOrbit`](@ref)s. Its parameters `θ` hold,
per orbit, the weight of each point and then the orbit's distinct coordinate values.
"""
struct BoxStructure
    D::Int
    orbits::Vector{BoxOrbit}
end
BoxStructure(mults::AbstractVector, D::Integer) = BoxStructure(Int(D), [BoxOrbit(m, D) for m in mults])
nunknowns(s::BoxStructure) = sum(nunknowns, s.orbits; init = 0)
npoints(s::BoxStructure) = sum(orbit_size, s.orbits; init = 0)
function param_offsets(s::BoxStructure)
    offs, k = Int[], 0
    for o in s.orbits
        push!(offs, k)
        k += nunknowns(o)
    end
    return offs
end
Base.show(io::IO, s::BoxStructure) = print(io, "BoxStructure(", s.D, ", ", [o.mult for o in s.orbits], ")")

# the representative's coordinates from the orbit's parameters
function representative(o::BoxOrbit, v::AbstractVector{S}) where {S}
    return [k == 0 ? zero(S) : v[k] for k in o.slots]
end

"Every signed permutation of `x`, without repeats."
function signed_images(x::AbstractVector)
    D = length(x)
    out = Vector{Vector{eltype(x)}}()
    seen = Set{Vector{eltype(x)}}()
    for σ in permutations_of(D), signs in 0:(2^D - 1)
        # a zero keeps its sign: -0.0 is a different key from 0.0 to `isequal`
        y = [((signs >> (j - 1)) & 1 == 1 && !iszero(x[σ[j]]) ? -x[σ[j]] : x[σ[j]]) for j in 1:D]
        y in seen && continue
        push!(seen, y)
        push!(out, y)
    end
    return out
end

"""
    expand(s::BoxStructure, θ) -> (nodes, weights)

The nodes and weights of the fully symmetric rule `(s, θ)` on `[-1, 1]^D`.
"""
function expand(s::BoxStructure, θ::AbstractVector{S}) where {S}
    nodes, ws = Vector{Vector{S}}(), S[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        imgs = signed_images(representative(o, θ[(off + 2):(off + nunknowns(o))]))
        length(imgs) == orbit_size(o) || throw(ArgumentError("orbit $(o.mult) degenerates: coinciding points"))
        append!(nodes, imgs)
        append!(ws, fill(θ[off + 1], length(imgs)))
    end
    return nodes, ws
end

"""
    box_margins(s, θ) -> (smallest weight, smallest distance of a coordinate from 0 and 1, smallest gap)

How far `(s, θ)` is from failing to be a positive rule with interior nodes in distinct
orbits: the third number is the smallest difference between two values of one orbit.
"""
function box_margins(s::BoxStructure, θ::AbstractVector{S}) where {S}
    wmin, cmin, gmin = typemax(S), typemax(S), typemax(S)
    for (o, off) in zip(s.orbits, param_offsets(s))
        wmin = min(wmin, θ[off + 1])
        v = θ[(off + 2):(off + nunknowns(o))]
        for i in eachindex(v)
            cmin = min(cmin, v[i], 1 - v[i])
            for j in (i + 1):length(v)
                gmin = min(gmin, abs(v[i] - v[j]))
            end
        end
    end
    return wmin, cmin, gmin
end

"`(s, θ)` with each orbit's values in a canonical order, and the orbits sorted."
function canonicalize(s::BoxStructure, θ::AbstractVector)
    blocks = Tuple{BoxOrbit,Vector{eltype(θ)}}[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        v = collect(θ[(off + 2):(off + nunknowns(o))])
        # values of equally repeated groups are interchangeable: largest first
        for m in unique(o.mult)
            idx = findall(==(m), o.mult)
            v[idx] = sort(abs.(v[idx]); rev = true)
        end
        push!(blocks, (o, vcat(θ[off + 1], abs.(v))))
    end
    sort!(blocks; by = b -> (orbit_size(b[1]), b[1].mult, Float64.(b[2][2:end])))
    return BoxStructure(s.D, [b[1] for b in blocks]), reduce(vcat, (b[2] for b in blocks); init = eltype(θ)[])
end

# --- the moment system ---------------------------------------------------------------------

"The exponent tuples `α₁ ≥ … ≥ α_D ≥ 0`, all even, of total degree at most `n`, by degree."
function box_exponents(D::Integer, n::Integer)
    out = Vector{Vector{Int}}()
    function rec(prefix, maxe, left)
        if length(prefix) == D
            push!(out, prefix)
            return
        end
        for e in 0:2:min(maxe, left)
            rec(vcat(prefix, e), e, left - e)
        end
    end
    rec(Int[], n, n)
    return sort!(out; by = α -> (sum(α), α))
end

"The distinct rearrangements of `α`."
distinct_permutations(α::AbstractVector) = unique([α[σ] for σ in permutations_of(length(α))])

"""
    BoxMomentSystem(structure, n, S)

The fully symmetric moment system for exactness to degree `n` on `[-1, 1]^D`, evaluated in
number type `S`: one equation per symmetrised product of even-degree orthonormal Legendre
polynomials. Callable: `sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)`.
"""
struct BoxMomentSystem{S}
    structure::BoxStructure
    n::Int
    exps::Vector{Vector{Int}}
    perms::Vector{Vector{Vector{Int}}}
    coef::Vector{S}
end
function BoxMomentSystem(structure::BoxStructure, n::Integer, ::Type{S}) where {S}
    exps = box_exponents(structure.D, n)
    perms = [distinct_permutations(α) for α in exps]
    return BoxMomentSystem{S}(structure, Int(n), exps, perms, [one(S) / sqrt(S(length(p))) for p in perms])
end
n_equations(sys::BoxMomentSystem) = length(sys.exps)
n_unknowns(sys::BoxMomentSystem) = nunknowns(sys.structure)

# orthonormal Legendre p_j(x) = √((2j+1)/2) P_j(x) and its derivative, j = 0 … n
function legendre_orthonormal(n::Integer, x::S) where {S}
    P, dP = zeros(S, n + 1), zeros(S, n + 1)
    P[1] = one(S)
    n >= 1 && (P[2] = x; dP[2] = one(S))
    for j in 1:(n - 1)
        P[j + 2] = ((2j + 1) * x * P[j + 1] - j * P[j]) / (j + 1)
        dP[j + 2] = dP[j] + (2j + 1) * P[j + 1]
    end
    c = [sqrt(S(2j + 1) / 2) for j in 0:n]
    return P .* c, dP .* c
end

function (sys::BoxMomentSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    s = sys.structure
    D, L = s.D, length(sys.exps)
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, nunknowns(s)) : zeros(S, 0, 0)
    for (o, off) in zip(s.orbits, param_offsets(s))
        w = S(θ[off + 1])
        sz = orbit_size(o)
        x = representative(o, S[θ[off + 1 + i] for i in 1:nparams(o)])
        tab = [legendre_orthonormal(sys.n, xk) for xk in x]
        for (k, α) in enumerate(sys.exps)
            φ = zero(S)
            g = zeros(S, D)                                  # ∂φ/∂x_j
            for σ in sys.perms[k]
                t = one(S)
                for j in 1:D
                    t *= tab[j][1][σ[j] + 1]
                end
                φ += t
                jacobian || continue
                for j in 1:D
                    u = tab[j][2][σ[j] + 1]
                    for l in 1:D
                        l == j || (u *= tab[l][1][σ[l] + 1])
                    end
                    g[j] += u
                end
            end
            φ *= sys.coef[k]
            r[k] += w * sz * φ
            jacobian || continue
            J[k, off + 1] += sz * φ
            for j in 1:D
                o.slots[j] == 0 && continue
                J[k, off + 1 + o.slots[j]] += w * sz * sys.coef[k] * g[j]
            end
        end
    end
    r[1] -= sqrt(S(2))^D                                     # ∫ φ_0 over the box; the rest integrate to 0
    return r, J
end
