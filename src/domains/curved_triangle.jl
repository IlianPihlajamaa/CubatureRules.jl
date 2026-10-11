# Curved boundary elements: the weakly singular kernel on a quadratic (six-node) triangle
# (notes/singular-bem.md, stage 6).
#
#     ∫_Γ f(y)/|y − x₀| dS_y = ∫_T̂ f(χ(ξ)) J(ξ)/|χ(ξ) − χ(ξ₀)| dξ,    Γ = χ(T̂), x₀ = χ(ξ₀),
#
# taken on the reference triangle T̂, where isoparametric shape functions are polynomials: a
# rule for the weight J(ξ)/|χ(ξ) − χ(ξ₀)| on T̂ is exact for every polynomial in ξ up to its
# degree.

"""
    QuadraticTriangle(p₁, p₂, p₃, p₁₂, p₂₃, p₃₁)

A curved triangle, the image of the reference triangle `Simplex{2}()` under the quadratic
Lagrange map `χ(ξ) = Σ Nᵢ(ξ) pᵢ` through its vertices `p₁, p₂, p₃` (at `ξ = (0, 0), (1, 0),
(0, 1)`) and the points `p₁₂, p₂₃, p₃₁` at the midpoints of its edges in `ξ`: the six-node
boundary element of isoparametric methods. Callable, `T(ξ) = χ(ξ)`; [`CurvedInverseDistance`](@ref)
puts the weakly singular kernel on it.
"""
struct QuadraticTriangle{D,T}
    p::SVector{6,SVector{D,T}}
end
function QuadraticTriangle(ps::Vararg{Any,6})
    D = length(ps[1])
    all(p -> length(p) == D, ps) || throw(DimensionMismatch("the six points need the same dimension"))
    T = promote_type((eltype(collect(p)) for p in ps)...)
    return QuadraticTriangle{D,T}(SVector{6,SVector{D,T}}(ntuple(i -> SVector{D,T}(Tuple(ps[i])), 6)))
end
Base.show(io::IO, t::QuadraticTriangle) = print(io, "QuadraticTriangle(", join((Tuple(p) for p in t.p), ", "), ")")

# barycentric coordinates and their (constant) gradients in ξ
_l(ξ) = (1 - ξ[1] - ξ[2], ξ[1], ξ[2])
# the arithmetic of a point of the element: the promotion of ξ's and the nodes' types, exact
# rationals widened to BigInt
_ctype(ξ, ::QuadraticTriangle{D,T}) where {D,T} = (S = promote_type(map(typeof, Tuple(ξ))..., T); S <: Union{Integer,Rational} ? Rational{BigInt} : S)
const _DL = ((-1, -1), (1, 0), (0, 1))

"`χ(ξ)`, in the arithmetic of `ξ` and the nodes promoted (exact rationals in `BigInt`)."
function (t::QuadraticTriangle{D})(ξ) where {D}
    S = _ctype(ξ, t)
    l1, l2, l3 = _l(S.(collect(ξ)))
    N = (l1 * (2l1 - 1), l2 * (2l2 - 1), l3 * (2l3 - 1), 4l1 * l2, 4l2 * l3, 4l3 * l1)
    return sum(N[i] * SVector{D,S}(t.p[i]) for i in 1:6)
end

"The `D × 2` Jacobian matrix `Dχ(ξ)`."
function jacobian_matrix(t::QuadraticTriangle{D}, ξ) where {D}
    S = _ctype(ξ, t)
    l = _l(S.(collect(ξ)))
    g = [SVector{2,S}(d) for d in _DL]
    dN = ((4l[1] - 1) * g[1], (4l[2] - 1) * g[2], (4l[3] - 1) * g[3],
          4(l[1] * g[2] + l[2] * g[1]), 4(l[2] * g[3] + l[3] * g[2]), 4(l[3] * g[1] + l[1] * g[3]))
    return sum(SVector{D,S}(t.p[i]) * transpose(dN[i]) for i in 1:6)
end

"The area element `J(ξ) = |∂₁χ × ∂₂χ|` (in two dimensions, `|det Dχ|`)."
function area_element(t::QuadraticTriangle{D}, ξ) where {D}
    A = jacobian_matrix(t, ξ)
    D == 2 && return abs(det(A))
    return norm(cross(A[:, 1], A[:, 2]))
end

"""
    second_term(t, v) = ½ vᵀ D²χ v

constant for a quadratic map, so that `χ(ξ₀ + s v) − χ(ξ₀) = s Dχ(ξ₀) v + s² second_term(t, v)`
exactly, without the cancellation of the difference.
"""
function second_term(t::QuadraticTriangle{D}, v::AbstractVector) where {D}
    S = _ctype(v, t)
    gv = [sum(d[k] * S(v[k]) for k in 1:2) for d in _DL]
    p = [SVector{D,S}(q) for q in t.p]
    return 2 * (gv[1]^2 * p[1] + gv[2]^2 * p[2] + gv[3]^2 * p[3]) +
           4 * (gv[1] * gv[2] * p[4] + gv[2] * gv[3] * p[5] + gv[3] * gv[1] * p[6])
end

"""
    CurvedInverseDistance(element, ξ₀)

The weakly singular kernel `1/|y − x₀|` on the curved triangle `element` (a
[`QuadraticTriangle`](@ref)), with `x₀ = element(ξ₀)` on it, as a weight on the reference
triangle: the weight of `WeightedDomain(Simplex{2}(), CurvedInverseDistance(element, ξ₀))` is

    w(ξ) = J(ξ) / |χ(ξ) − χ(ξ₀)|,

`J` the area element, so that `∫ g(ξ) w(ξ) dξ = ∫_Γ g(χ⁻¹(y))/|y − x₀| dS_y`. A rule for it
is exact when `g` is a polynomial in `ξ` of its degree — the shape functions of isoparametric
elements — and its nodes are in `ξ`: the surface points are `element.(nodes(r))`.

```julia
Γ = QuadraticTriangle((1, 0, 0), (0, 1, 0), (0, 0, 1), (0.7071, 0.7071, 0), (0, 0.7071, 0.7071), (0.7071, 0, 0.7071))
dom = WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ, (1//3, 1//3)))
r = rule(dom; degree = 7)
integrate(ξ -> 1 - ξ[1] - ξ[2], r)          # ∫_Γ N₁/|y − x₀| dS, N₁ the linear shape function
```

`ξ₀` may be a vertex, a point of an edge or an interior point.

    CurvedInverseDistance(element; x0)

The same kernel for a point `x₀` given in space, on the surface or off it: the weight
`J(ξ)/|χ(ξ) − x₀|`. Off the surface the integral is near-singular as `x₀` comes close; the
rule is built about the point of the element nearest `x₀`, found by Newton's method on
`|χ(ξ) − x₀|²` (inside the element, along its edges and at its vertices). A point within
rounding of the surface is taken as on it.

Rules come from [`DuffyCurved`](@ref).
"""
struct CurvedInverseDistance{E<:QuadraticTriangle,P<:Real,X}
    element::E
    ξ0::SVector{2,P}            # the parameter point of x₀ when that is given
    x0::X                       # x₀ in space when that is given, else nothing
end
CurvedInverseDistance(e::QuadraticTriangle, ξ0::Union{Tuple,AbstractVector}) =
    CurvedInverseDistance(e, SVector{2}(promote(ξ0...)), nothing)
function CurvedInverseDistance(e::QuadraticTriangle{D}; x0) where {D}
    length(x0) == D || throw(DimensionMismatch("x₀ needs $D coordinates, like the element"))
    all(isfinite, x0) || throw(ArgumentError("the point must be finite, got $(Tuple(x0))"))
    return CurvedInverseDistance(e, SVector{2,Int}(0, 0), SVector{D}(promote(x0...)))
end
Base.show(io::IO, w::CurvedInverseDistance) =
    w.x0 === nothing ? print(io, "CurvedInverseDistance(", w.element, ", ", Tuple(w.ξ0), ")") :
                       print(io, "CurvedInverseDistance(", w.element, "; x0 = ", Tuple(w.x0), ")")

"The reference triangle carrying a [`CurvedInverseDistance`](@ref) kernel."
const CurvedKernelDomain = WeightedDomain{2,<:Any,<:Simplex{2},<:CurvedInverseDistance}

"""
    nearest_point(element, X) -> (ξ, distance)

The point of the reference triangle whose image is nearest `X`, in the arithmetic of `X`:
Newton's method on `½|χ(ξ) − X|²` from the centroid and from halfway between it and each vertex
(with the constant second derivatives of the quadratic map), Newton along each edge,
and the vertices; the closest of the candidates inside the closed triangle.
"""
function nearest_point(el::QuadraticTriangle{D}, X::SVector{D,S}) where {D,S}
    F(ξ) = (d = el(ξ) - X; dot(d, d))
    # the Hessians of the components, from ½ vᵀ D²χ v
    e1, e2 = SVector{2,S}(1, 0), SVector{2,S}(0, 1)
    Q1, Q2, Q12 = second_term(el, e1), second_term(el, e2), second_term(el, e1 + e2)
    H11, H22 = 2Q1, 2Q2
    H12 = Q12 - Q1 - Q2
    tol = ldexp(one(S), -(precision(S) - 8))
    inside(ξ) = all(isfinite, ξ) && ξ[1] >= 0 && ξ[2] >= 0 && ξ[1] + ξ[2] <= 1
    function newton2(ξ)
        for _ in 1:80
            d = el(ξ) - X
            A = jacobian_matrix(el, ξ)
            g = transpose(A) * d
            M = transpose(A) * A + SMatrix{2,2,S,4}(dot(d, H11), dot(d, H12), dot(d, H12), dot(d, H22))
            Δ = M \ g
            ξ = ξ - Δ
            (all(isfinite, ξ) && norm(ξ) < 10) || return nothing
            norm(Δ) <= tol && return ξ
        end
        return ξ
    end
    function newton1(p, q)
        u = q - p
        τ = one(S) / 2
        for _ in 1:80
            ξ = p + τ * u
            d = el(ξ) - X
            Au = jacobian_matrix(el, ξ) * u
            h = dot(Au, Au) + 2dot(d, second_term(el, u))       # second derivative in τ
            Δ = h > 0 ? dot(Au, d) / h : (dot(Au, d) > 0 ? τ : τ - 1)
            τ = clamp(τ - Δ, zero(S), one(S))
            abs(Δ) <= tol && break
        end
        return p + τ * u
    end
    V = (SVector{2,S}(0, 0), SVector{2,S}(1, 0), SVector{2,S}(0, 1))
    cands = collect(V)
    for (p, q) in ((V[1], V[2]), (V[2], V[3]), (V[3], V[1]))
        push!(cands, newton1(p, q))
    end
    for start in (SVector{2,S}(2, 2), SVector{2,S}(1, 1), SVector{2,S}(4, 1), SVector{2,S}(1, 4)) ./ 6
        ξ = newton2(start)
        ξ !== nothing && inside(ξ) && push!(cands, ξ)
    end
    best = cands[argmin(F.(cands))]
    return best, sqrt(F(best))
end

"""
    curved_geometry(dom; bits = 256) -> (; place, x0, patches, moved, X, H)

How the reference triangle is cut, as for the flat case (see `kernel_geometry`): at `ξ₀` when
that is given, and at the parameter point of the element nearest `x₀` when `x₀` is given in
space, computed at `bits`. `place` is `:vertex`, `:edge` or `:interior` for `x₀` on the surface
(within rounding: `moved` says how far), and `:near` off it; `X` is then `x₀` and `H` its
distance from the element, as `BigFloat`s.
"""
function curved_geometry(dom::CurvedKernelDomain; bits::Integer = 256)
    w = dom.weight
    if w.x0 === nothing
        g = kernel_geometry(WeightedDomain(dom.base, InverseDistance(w.ξ0)))
        return (; place = g.place, x0 = g.x0, patches = g.patches, moved = g.moved, X = nothing, H = big(0.0))
    end
    D = length(w.x0)
    ξc, H, X = with_bits(bits) do
        X = SVector{D,BigFloat}(map(BigFloat, w.x0))
        ξ, h = nearest_point(w.element, X)
        ξ, h, X
    end
    g = kernel_geometry(WeightedDomain(dom.base, InverseDistance(Tuple(ξc))))
    scale = with_bits(() -> max(one(BigFloat), maximum(abs, X)), bits)
    tol = 64 * maximum(x -> BigFloat(value_ulp(x)), collect(w.x0)) + ldexp(scale, -(bits - 32))
    H <= tol && return (; place = g.place, x0 = g.x0, patches = g.patches, moved = Float64(H), X = nothing, H = big(0.0))
    return (; place = :near, x0 = g.x0, patches = g.patches, moved = 0.0, X, H)
end

"""
`∫_T̂ J(ξ)/|χ(ξ) − χ(ξ₀)| dξ`, the area-weighted `∫_Γ dS/|y − x₀|`, in the ambient `BigFloat`
precision (see `polar_curved_integrals`).
"""
measure(d::CurvedKernelDomain) = only(polar_curved_integrals(d, 1, BigFloat, (out, ξ) -> (out[1] = one(BigFloat))))
