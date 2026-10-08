# Fully symmetric rules on the pyramid (PLAN §6 Tier 3): the reference pyramid, base [-1, 1]²
# at z = 0 and apex (0, 0, 1), is invariant under the 8 symmetries of its square base acting on
# (x, y), with z left alone (the group C₄ᵥ).
#
# An orbit is a square orbit (symmetry/box.jl) at some height. It is written in the collapsed
# coordinates ξ = x/(1 − z), η = y/(1 − z), in which every cross-section is [-1, 1]², so an
# orbit's parameters are its distinct |ξ| values, each in (0, 1), and its height z in (0, 1).
#
# The moment equations use the collapsed-coordinate basis
#
#     φ = S_ij(ξ, η) (1 − z)^(i+j) q_k(z),   S_ij = (p_i(ξ) p_j(η) + p_j(ξ) p_i(η)) / √2 (or p_i p_i),
#
# with p orthonormal Legendre, i ≥ j both even, and q_k orthonormal on [0, 1] for the weight
# (1 − z)^(2(i+j)+2). p_i(ξ) (1 − z)^i is a polynomial in x and z, so φ is a polynomial of
# degree i + j + k; with dx dy = (1 − z)² dξ dη these are orthonormal on the pyramid, and the
# symmetrised products are exactly the invariant ones. Only the constant has a nonzero integral.

"""
    PyramidOrbit(base)

An orbit type of the pyramid's symmetry group: the square orbit `base` (a [`BoxOrbit`](@ref)
with `D = 2`) in the collapsed coordinates, at a height that is one more parameter.
"""
struct PyramidOrbit
    base::BoxOrbit
end
PyramidOrbit(mult::AbstractVector) = PyramidOrbit(BoxOrbit(mult, 2))
Base.:(==)(a::PyramidOrbit, b::PyramidOrbit) = a.base == b.base
Base.hash(o::PyramidOrbit, h::UInt) = hash(o.base, hash(:PyramidOrbit, h))
nparams(o::PyramidOrbit) = nparams(o.base) + 1                 # the base's values, then z
nunknowns(o::PyramidOrbit) = 1 + nparams(o)
orbit_size(o::PyramidOrbit) = orbit_size(o.base)
pyramid_orbit_types() = [PyramidOrbit(o) for o in box_orbit_types(2)]

"""
    PyramidStructure(orbits)

A fully symmetric rule on the reference pyramid as a list of [`PyramidOrbit`](@ref)s. Its
parameters hold, per orbit, the weight of each point, the distinct `|ξ|` values and the height.
"""
struct PyramidStructure
    orbits::Vector{PyramidOrbit}
end
PyramidStructure(mults::AbstractVector{<:AbstractVector}) = PyramidStructure([PyramidOrbit(m) for m in mults])
nunknowns(s::PyramidStructure) = sum(nunknowns, s.orbits; init = 0)
npoints(s::PyramidStructure) = sum(orbit_size, s.orbits; init = 0)
function param_offsets(s::PyramidStructure)
    offs, k = Int[], 0
    for o in s.orbits
        push!(offs, k)
        k += nunknowns(o)
    end
    return offs
end
Base.show(io::IO, s::PyramidStructure) = print(io, "PyramidStructure(", [o.base.mult for o in s.orbits], ")")

"""
    expand(s::PyramidStructure, θ) -> (nodes, weights)

The nodes (Cartesian, on the reference pyramid) and weights of the rule `(s, θ)`.
"""
function expand(s::PyramidStructure, θ::AbstractVector{S}) where {S}
    nodes, ws = Vector{Vector{S}}(), S[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        k = nparams(o.base)
        z = θ[off + 2 + k]
        imgs = signed_images(representative(o.base, θ[(off + 2):(off + 1 + k)]))
        length(imgs) == orbit_size(o) || throw(ArgumentError("orbit $(o.base.mult) degenerates: coinciding points"))
        for ξ in imgs
            push!(nodes, [ξ[1] * (1 - z), ξ[2] * (1 - z), z])
            push!(ws, θ[off + 1])
        end
    end
    return nodes, ws
end

"The invariant basis terms `(i, j, k)` of degree at most `n`, by degree."
pyramid_terms(n::Integer) = sort!([(i, j, k) for i in 0:2:n for j in 0:2:i for k in 0:(n - i - j)];
                                  by = t -> (sum(t), t[1], t[3]))

"""
    PyramidMomentSystem(structure, n, S)

The fully symmetric moment system for exactness to degree `n` on the reference pyramid, in
the orthonormal invariant polynomials described at the top of this file. Callable:
`sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)`.
"""
struct PyramidMomentSystem{S}
    structure::PyramidStructure
    n::Int
    terms::Vector{NTuple{3,Int}}
end
PyramidMomentSystem(s::PyramidStructure, n::Integer, ::Type{S}) where {S} =
    PyramidMomentSystem{S}(s, Int(n), pyramid_terms(n))
n_equations(sys::PyramidMomentSystem) = length(sys.terms)
n_unknowns(sys::PyramidMomentSystem) = nunknowns(sys.structure)

function (sys::PyramidMomentSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    s, n = sys.structure, sys.n
    L = length(sys.terms)
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, nunknowns(s)) : zeros(S, 0, 0)
    for (o, off) in zip(s.orbits, param_offsets(s))
        w = S(θ[off + 1])
        sz = orbit_size(o)
        k = nparams(o.base)
        v = S[θ[off + 1 + i] for i in 1:k]
        z = S(θ[off + 2 + k])
        ξ = representative(o.base, v)
        pξ, dpξ = legendre_orthonormal(n, ξ[1])
        pη, dpη = legendre_orthonormal(n, ξ[2])
        t = 2z - 1
        for (row, (i, j, kk)) in enumerate(sys.terms)
            a = 2(i + j) + 2
            q, dq = jacobi_orthonormal_d(kk, a, 0, t)
            c = sqrt(ldexp(one(S), a + 1))                     # q_k(z) = 2^((a+1)/2) p_k(2z − 1)
            g = (1 - z)^(i + j) * c * q[kk + 1]
            dg = (i + j == 0 ? zero(S) : -(i + j) * (1 - z)^(i + j - 1) * c * q[kk + 1]) + (1 - z)^(i + j) * c * 2 * dq[kk + 1]
            # the symmetrised product and its derivatives in ξ and η
            if i == j
                Sij, dSξ, dSη = pξ[i + 1] * pη[i + 1], dpξ[i + 1] * pη[i + 1], pξ[i + 1] * dpη[i + 1]
            else
                h = one(S) / sqrt(S(2))
                Sij = h * (pξ[i + 1] * pη[j + 1] + pξ[j + 1] * pη[i + 1])
                dSξ = h * (dpξ[i + 1] * pη[j + 1] + dpξ[j + 1] * pη[i + 1])
                dSη = h * (pξ[i + 1] * dpη[j + 1] + pξ[j + 1] * dpη[i + 1])
            end
            φ = Sij * g
            r[row] += w * sz * φ
            jacobian || continue
            J[row, off + 1] += sz * φ
            for q_ in 1:k                                      # the base's values: which of ξ, η carry them
                d = (o.base.slots[1] == q_ ? dSξ : zero(S)) + (o.base.slots[2] == q_ ? dSη : zero(S))
                J[row, off + 1 + q_] += w * sz * d * g
            end
            J[row, off + 2 + k] += w * sz * Sij * dg
        end
    end
    r[1] -= 2 / sqrt(S(3))                                     # ∫ φ₀₀₀ = (√3/2) · 4/3
    return r, J
end

"`(smallest weight, smallest distance of a parameter from the edge of its range, smallest gap)`"
function pyramid_margins(s::PyramidStructure, θ::AbstractVector{S}) where {S}
    wmin, cmin, gmin = typemax(S), typemax(S), typemax(S)
    for (o, off) in zip(s.orbits, param_offsets(s))
        wmin = min(wmin, θ[off + 1])
        k = nparams(o.base)
        v = θ[(off + 2):(off + 1 + k)]
        z = θ[off + 2 + k]
        cmin = min(cmin, z, 1 - z)
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
function canonicalize(s::PyramidStructure, θ::AbstractVector)
    blocks = Tuple{PyramidOrbit,Vector{eltype(θ)}}[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        k = nparams(o.base)
        v = abs.(collect(θ[(off + 2):(off + 1 + k)]))
        for m in unique(o.base.mult)                          # equally repeated values are interchangeable
            idx = findall(==(m), o.base.mult)
            v[idx] = sort(v[idx]; rev = true)
        end
        push!(blocks, (o, vcat(θ[off + 1], v, θ[off + 2 + k])))
    end
    sort!(blocks; by = b -> (orbit_size(b[1]), b[1].base.mult, Float64.(b[2][2:end])))
    return PyramidStructure([b[1] for b in blocks]), reduce(vcat, (b[2] for b in blocks); init = eltype(θ)[])
end
