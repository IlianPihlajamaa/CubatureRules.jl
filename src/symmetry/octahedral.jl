# Orbit algebra for the octahedral group O_h on the sphere (PLAN §6 Tier 3, v0.4), the
# symmetry of the Lebedev rules. Where the simplex machinery has S_N permuting barycentric
# coordinates, this has all 48 signed permutations of Cartesian ones.
#
# The orbits are the classical six:
#
#     a₁   (1,0,0)                    6 points,  no parameter
#     a₂   (1,1,0)/√2                12 points,  no parameter
#     a₃   (1,1,1)/√3                 8 points,  no parameter
#     b    (l,l,m), m² = 1-2l²       24 points,  1 parameter
#     c    (p,q,0), q² = 1-p²        24 points,  1 parameter
#     d    (r,s,u), u² = 1-r²-s²     48 points,  2 parameters
#
# each carrying one weight besides.
#
# Two facts make the moment system small and rational, and both are worth stating plainly.
#
# First, a rule assembled from whole orbits is O_h-invariant, so its residual against any
# harmonic whose group average vanishes is *identically* zero — no information at all. At
# degree 131 that is 17424 harmonics of which a few hundred say anything. The system is
# therefore written in the invariants from the start. The invariant ring of O_h on R³ is
# free on p₂ = Σxᵢ², p₄ = Σxᵢ⁴ and p₆ = x²y²z², and p₂ = 1 on the sphere, so the invariant
# polynomials of degree ≤ n restricted to the sphere are spanned by p₄^a p₆^b with
# 4a + 6b ≤ n. (Checked against the group-averaged monomials in the test suite, not assumed.)
#
# Second, an invariant is *constant on an orbit*, so evaluating it costs one point rather
# than 48, and p₄ and p₆ are polynomials in the orbit parameters with the constraint
# substituted — no square roots appear in the equations at all, only in the nodes.

"""
    OctahedralOrbit(kind)

One `O_h`-orbit type on the sphere: `:a1`, `:a2`, `:a3`, `:b`, `:c` or `:d`, as in the table
at the top of this file.
"""
struct OctahedralOrbit
    kind::Symbol
    function OctahedralOrbit(kind::Symbol)
        kind in (:a1, :a2, :a3, :b, :c, :d) ||
            throw(ArgumentError("unknown octahedral orbit $(kind); expected one of :a1, :a2, :a3, :b, :c, :d"))
        return new(kind)
    end
end

Base.show(io::IO, o::OctahedralOrbit) = print(io, "OctahedralOrbit(:", o.kind, ")")

"Points in the orbit: 6, 12, 8, 24, 24 or 48."
orbit_size(o::OctahedralOrbit) = o.kind === :a1 ? 6 : o.kind === :a2 ? 12 : o.kind === :a3 ? 8 :
                                 o.kind === :d ? 48 : 24
"Free parameters of the orbit, the weight excluded."
nparams(o::OctahedralOrbit) = o.kind in (:a1, :a2, :a3) ? 0 : o.kind === :d ? 2 : 1
"Unknowns contributed to a moment system: the weight and the free parameters."
nunknowns(o::OctahedralOrbit) = 1 + nparams(o)

"""
    orbit_representative(o, p, S) -> SVector{3,S}

One point of the orbit, from its parameters `p`. Every `O_h`-invariant takes the same value
at every point of the orbit, so this is all the moment system ever needs.
"""
function orbit_representative(o::OctahedralOrbit, p, ::Type{S}) where {S}
    k = o.kind
    k === :a1 && return SVector{3,S}(one(S), zero(S), zero(S))
    k === :a2 && return SVector{3,S}(one(S) / sqrt(S(2)), one(S) / sqrt(S(2)), zero(S))
    k === :a3 && return (c = one(S) / sqrt(S(3)); SVector{3,S}(c, c, c))
    if k === :b
        l = S(p[1])
        m2 = one(S) - 2l^2
        m2 >= 0 || throw(DomainError(l, "orbit b needs 2l² ≤ 1"))
        return SVector{3,S}(l, l, sqrt(m2))
    elseif k === :c
        q = S(p[1])
        r2 = one(S) - q^2
        r2 >= 0 || throw(DomainError(q, "orbit c needs p² ≤ 1"))
        return SVector{3,S}(q, sqrt(r2), zero(S))
    else
        r, s = S(p[1]), S(p[2])
        u2 = one(S) - r^2 - s^2
        u2 >= 0 || throw(DomainError((r, s), "orbit d needs r² + s² ≤ 1"))
        return SVector{3,S}(r, s, sqrt(u2))
    end
end

"""
    signed_orbit(v) -> Vector{SVector{3}}

Every distinct image of `v` under the 48 signed permutations of the coordinates. Permuting
and negating are exact, so duplicates are bitwise identical and drop out — but a sign on a
zero coordinate is not, so zeros are never negated.
"""
function signed_orbit(v::SVector{3,S}) where {S}
    pts = Vector{SVector{3,S}}()
    for p in permutations_of(3), s1 in (1, -1), s2 in (1, -1), s3 in (1, -1)
        a, b, c = v[p[1]], v[p[2]], v[p[3]]
        w = SVector{3,S}(iszero(a) ? a : s1 * a, iszero(b) ? b : s2 * b, iszero(c) ? c : s3 * c)
        push!(pts, w)
    end
    return unique(pts)
end

orbit_points(o::OctahedralOrbit, p, ::Type{S}) where {S} = signed_orbit(orbit_representative(o, p, S))

"""
    OctahedralStructure(orbits)

A list of `O_h`-orbit types: the shape of a Lebedev-like rule, before any numbers.
"""
struct OctahedralStructure
    orbits::Vector{OctahedralOrbit}
end
OctahedralStructure(kinds::AbstractVector{Symbol}) = OctahedralStructure(OctahedralOrbit.(kinds))

npoints(s::OctahedralStructure) = sum(orbit_size, s.orbits; init = 0)
nunknowns(s::OctahedralStructure) = sum(nunknowns, s.orbits; init = 0)
Base.show(io::IO, s::OctahedralStructure) =
    print(io, "OctahedralStructure(", [o.kind for o in s.orbits], ", ", npoints(s), " points, ",
          nunknowns(s), " unknowns)")

"Offsets of each orbit's parameter block in θ = [w₁, p₁…, w₂, p₂…, …]."
function param_offsets(s::OctahedralStructure)
    offs = Vector{Int}(undef, length(s.orbits))
    k = 0
    for (i, o) in enumerate(s.orbits)
        offs[i] = k
        k += nunknowns(o)
    end
    return offs
end

"""
    expand(s::OctahedralStructure, θ) -> (nodes, weights)

The nodes and weights of the rule described by `θ = [w₁, p₁…, w₂, p₂…, …]`, where each `w`
is the weight *per point*.
"""
function expand(s::OctahedralStructure, θ::AbstractVector{S}) where {S}
    length(θ) == nunknowns(s) || throw(ArgumentError("expected $(nunknowns(s)) unknowns, got $(length(θ))"))
    offs = param_offsets(s)
    xs = Vector{SVector{3,S}}()
    ws = Vector{S}()
    for (i, o) in enumerate(s.orbits)
        w = θ[offs[i] + 1]
        p = θ[(offs[i] + 2):(offs[i] + nunknowns(o))]
        pts = orbit_points(o, p, S)
        length(pts) == orbit_size(o) ||
            throw(ArgumentError("orbit $(o.kind) degenerated to $(length(pts)) points at these parameters"))
        append!(xs, pts)
        append!(ws, fill(w, length(pts)))
    end
    return xs, ws
end

# --- the invariants ----------------------------------------------------------------------
#
# p₄ and p₆ at an orbit's representative, as functions of the parameters with the constraint
# already substituted, together with their derivatives. Polynomials throughout: the square
# roots live in the nodes, not in the equations.

"""
    orbit_invariants(o, p, S) -> (p₄, p₆, ∂p₄, ∂p₆)

The two invariants at the orbit's points and their gradients with respect to the orbit's
free parameters (empty for the parameterless orbits).
"""
function orbit_invariants(o::OctahedralOrbit, p, ::Type{S}) where {S}
    k = o.kind
    if k === :a1                       # (1,0,0)
        return one(S), zero(S), S[], S[]
    elseif k === :a2                   # (1,1,0)/√2 → p₄ = 2·(1/4)
        return one(S) / 2, zero(S), S[], S[]
    elseif k === :a3                   # (1,1,1)/√3 → p₄ = 3·(1/9), p₆ = (1/3)³
        return one(S) / 3, one(S) / 27, S[], S[]
    elseif k === :b                    # (l,l,m), m² = 1-2l²
        l = S(p[1])
        l2 = l^2
        m2 = one(S) - 2l2
        p4 = 2l2^2 + m2^2
        p6 = l2^2 * m2
        # d/dl: dl² = 2l, dm² = -4l
        dp4 = 4l2 * 2l + 2 * m2 * (-4l)
        dp6 = 2l2 * 2l * m2 + l2^2 * (-4l)
        return p4, p6, S[dp4], S[dp6]
    elseif k === :c                    # (p,q,0), q² = 1-p²
        a = S(p[1])
        a2 = a^2
        b2 = one(S) - a2
        p4 = a2^2 + b2^2
        dp4 = 4a2 * a + 2 * b2 * (-2a)
        return p4, zero(S), S[dp4], S[zero(S)]
    else                               # (r,s,u), u² = 1-r²-s²
        r, s = S(p[1]), S(p[2])
        r2, s2 = r^2, s^2
        u2 = one(S) - r2 - s2
        p4 = r2^2 + s2^2 + u2^2
        p6 = r2 * s2 * u2
        dp4_r = 4r2 * r + 2 * u2 * (-2r)
        dp4_s = 4s2 * s + 2 * u2 * (-2s)
        dp6_r = 2r * s2 * u2 + r2 * s2 * (-2r)
        dp6_s = r2 * 2s * u2 + r2 * s2 * (-2s)
        return p4, p6, S[dp4_r, dp4_s], S[dp6_r, dp6_s]
    end
end

"""
    invariant_exponents(n) -> Vector{Tuple{Int,Int}}

The exponents `(a, b)` of the invariant basis `p₄^a p₆^b` with `4a + 6b ≤ n`: a basis of the
`O_h`-invariant polynomials of degree `≤ n` restricted to the sphere. The constant `(0,0)`
comes first.
"""
invariant_exponents(n::Integer) =
    sort!([(a, b) for a in 0:(n ÷ 4) for b in 0:(n ÷ 6) if 4a + 6b <= n]; by = t -> (4t[1] + 6t[2], t[2]))

"""
    invariant_moment(a, b) -> BigFloat

`∫_{S²} p₄^a p₆^b dσ`, by multinomial expansion into monomials whose moments are known
exactly:  `p₄^a = Σ (a choose i,j,k) x^{4i} y^{4j} z^{4k}`, and `p₆^b` shifts every exponent
by `2b`.
"""
function invariant_moment(a::Integer, b::Integer)
    total = zero(BigFloat)
    for i in 0:a, j in 0:(a - i)
        k = a - i - j
        coeff = BigFloat(factorial(big(a)) ÷ (factorial(big(i)) * factorial(big(j)) * factorial(big(k))))
        total += coeff * sphere_moment(3, (4i + 2b, 4j + 2b, 4k + 2b))
    end
    return total
end

# --- the moment system -------------------------------------------------------------------

"""
    OctahedralMomentSystem(structure, n, S)

The reduced moment system for exactness to degree `n` of an `O_h`-symmetric rule on the
sphere: one equation per invariant `p₄^a p₆^b` with `4a + 6b ≤ n`, in number type `S`.
Callable, `sys(θ) -> (r, J)`, which is what [`gauss_newton`](@ref) consumes.

The equations are

    Σ_orbits wₒ nₒ p₄(θₒ)^a p₆(θₒ)^b  =  ∫_{S²} p₄^a p₆^b dσ,

with `nₒ` the orbit size — an invariant is constant on an orbit, so each orbit contributes
through one representative rather than 48 points.
"""
struct OctahedralMomentSystem{S}
    structure::OctahedralStructure
    n::Int
    exps::Vector{Tuple{Int,Int}}
    rhs::Vector{S}
end

function OctahedralMomentSystem(structure::OctahedralStructure, n::Integer, ::Type{S}) where {S}
    exps = invariant_exponents(n)
    rhs = [S(invariant_moment(a, b)) for (a, b) in exps]
    return OctahedralMomentSystem{S}(structure, Int(n), exps, rhs)
end

n_equations(sys::OctahedralMomentSystem) = length(sys.exps)
n_unknowns(sys::OctahedralMomentSystem) = nunknowns(sys.structure)

function (sys::OctahedralMomentSystem{S})(θ::AbstractVector) where {S}
    L, p = n_equations(sys), n_unknowns(sys)
    length(θ) == p || throw(ArgumentError("expected $p unknowns, got $(length(θ))"))
    r = -copy(sys.rhs)
    J = zeros(S, L, p)
    offs = param_offsets(sys.structure)
    for (i, o) in enumerate(sys.structure.orbits)
        w = S(θ[offs[i] + 1])
        par = θ[(offs[i] + 2):(offs[i] + nunknowns(o))]
        n_o = S(orbit_size(o))
        p4, p6, dp4, dp6 = orbit_invariants(o, par, S)
        for (k, (a, b)) in enumerate(sys.exps)
            ψ = p4^a * p6^b
            r[k] += w * n_o * ψ
            J[k, offs[i] + 1] = n_o * ψ
            for j in 1:nparams(o)
                # ∂(p₄^a p₆^b), with the b = 0 and a = 0 cases dropped rather than divided
                d = zero(S)
                a > 0 && (d += a * p4^(a - 1) * p6^b * dp4[j])
                b > 0 && (d += b * p4^a * p6^(b - 1) * dp6[j])
                J[k, offs[i] + 1 + j] = w * n_o * d
            end
        end
    end
    return r, J
end

"""
    refine_octahedral(structure, n, θ64, bits; cancel) -> (θ, result, guard_bits)

Refine a `Float64` seed of an `O_h`-symmetric sphere rule to `bits` bits, the same way
[`refine_symmetric`](@ref) does on a simplex: guard digits from the condition number
measured at the seed, Gauss–Newton at `bits + guard`, and one re-run if the condition
number met along the way asks for more.
"""
function refine_octahedral(structure::OctahedralStructure, n::Integer, θ64::Vector{Float64},
                           bits::Integer; cancel = nothing)
    sys64 = OctahedralMomentSystem(structure, n, Float64)
    r64, J64 = sys64(θ64)
    κ0 = lsq_step(J64, r64; rank_rtol = 1e-14)[2]
    guard = guard_bits_from_cond(κ0)
    for attempt in 1:2
        wbits = bits + guard
        res = with_bits(wbits) do
            sys = OctahedralMomentSystem(structure, n, BigFloat)
            gauss_newton(sys, BigFloat.(θ64); step_tol = ldexp(BigFloat(1), -(bits + 16)),
                         res_floor = ldexp(BigFloat(1), -(wbits - 12)),
                         rank_rtol = ldexp(BigFloat(1), -(wbits ÷ 2)), maxiter = 60, cancel)
        end
        needed = guard_bits_from_cond(res.cond_max)
        if needed <= guard || attempt == 2
            return res.θ, res, guard
        end
        guard = needed
    end
end
