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

`ξ₀` may be a vertex, a point of an edge or an interior point; rules come from
[`DuffyCurved`](@ref).
"""
struct CurvedInverseDistance{E<:QuadraticTriangle,P<:Real}
    element::E
    ξ0::SVector{2,P}
end
CurvedInverseDistance(e::QuadraticTriangle, ξ0::Union{Tuple,AbstractVector}) =
    CurvedInverseDistance(e, SVector{2}(promote(ξ0...)))
Base.show(io::IO, w::CurvedInverseDistance) = print(io, "CurvedInverseDistance(", w.element, ", ", Tuple(w.ξ0), ")")

"The reference triangle carrying a [`CurvedInverseDistance`](@ref) kernel."
const CurvedKernelDomain = WeightedDomain{2,<:Any,<:Simplex{2},<:CurvedInverseDistance}

"How the reference triangle is cut at `ξ₀`: the geometry of the flat case (see `kernel_geometry`)."
curved_geometry(dom::CurvedKernelDomain) = kernel_geometry(WeightedDomain(dom.base, InverseDistance(dom.weight.ξ0)))

"""
`∫_T̂ J(ξ)/|χ(ξ) − χ(ξ₀)| dξ`, the area-weighted `∫_Γ dS/|y − x₀|`, in the ambient `BigFloat`
precision (see `polar_curved_integrals`).
"""
measure(d::CurvedKernelDomain) = only(polar_curved_integrals(d, 1, BigFloat, (out, ξ) -> (out[1] = one(BigFloat))))
