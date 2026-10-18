# Transporting and combining rules, each with a stated effect on the claim (PLAN §2.3):
#
#     map_to(r, dom)        affine only          claim preserved
#     subdivide(r, dom, n)  affine pieces        claim preserved (per cell, hence globally)
#     transform(r, φ, Jφ)   general φ            → NoClaim, unless the caller asserts one
#     duffy(r)              radial grading       → NoClaim

"Run `f` at the precision of the rule's BigFloats (so no conversion reads ambient precision)."
_at_rule_precision(f, r::QuadratureRule) = numbertype(r) === BigFloat ? with_bits(f, _rule_precision(r)) : f()

"""
    map_to(r, domain)

The affine image of rule `r` on `domain` — a domain of the same kind and dimension, or a
triangle in space for a rule on a triangle in the plane — with the weights scaled by the
Jacobian. The claim is preserved. There is deliberately no way to pass a nonlinear map here:
use [`transform`](@ref), which reports `NoClaim`.
"""
function map_to(r::QuadratureRule, dom::Domain)
    S = _maptype(numbertype(r), _coordtype(dom))
    kind = dom isa Union{Sphere,Ball} ? "similarity" : "affine"
    return _at_rule_precision(r) do
        A, b, J = _affine(r.domain, dom, S)
        xs = [A * _mapnode(S, x) + b for x in r.nodes]
        ws = [S(w) * J for w in r.weights]
        QuadratureRule(xs, ws, convert_domain(S, dom), r.exactness,
                       with_step(r.provenance, "map_to: $kind image onto $(dom)"), r.certificate)
    end
end

# The pairs of domains map_to and integrate over a domain transport between.
_mappable(::Domain, ::Domain) = false
_mappable(::Interval, ::Interval) = true
_mappable(::Simplex{D}, ::Simplex{D}) where {D} = true
_mappable(::Simplex{2}, ::SurfaceTriangle) = true
_mappable(::Orthotope{D}, ::Orthotope{D}) where {D} = true
_mappable(::Wedge, ::Wedge) = true
_mappable(::Pyramid, ::Pyramid) = true
_mappable(::Ball{D}, ::Ball{D}) where {D} = true
_mappable(::Sphere{D}, ::Sphere{D}) where {D} = true
_mappable(a::JacobiInterval, b::JacobiInterval) = a.weight.α == b.weight.α && a.weight.β == b.weight.β

_mapnode(::Type{S}, x::Number) where {S} = S(x)
_mapnode(::Type{S}, x::SVector{D}) where {S,D} = SVector{D,S}(map(S, x))

# (A, b, J) of the affine map from `src` onto `dst`, in type S: frame(dst) ∘ frame(src)⁻¹
# (see `affine_frame`), and the factor by which it scales the measure. A reference `src`
# needs no inversion; allocation-free for isbits S.
@inline function _affine(src::Domain, dst::Domain, ::Type{S}) where {S}
    _mappable(src, dst) ||
        throw(ArgumentError("map_to transports a rule affinely between domains of the same kind; " *
                            "cannot map a rule on $(src) onto $(dst)"))
    Ad, bd, Jd = affine_frame(dst, S)
    isreference(src) && return Ad, bd, Jd
    As, bs, Js = affine_frame(src, S)
    A = Ad / As
    return A, bd - A * bs, _composed_jacobian(dst, A, Jd, Js)
end
_composed_jacobian(::Domain, A, Jd, Js) = abs(det(A))
_composed_jacobian(::Interval, A, Jd, Js) = abs(A)
_composed_jacobian(::Union{Sphere,SurfaceTriangle}, A, Jd, Js) = Jd / Js    # a surface: the ratio of area elements
_composed_jacobian(::JacobiInterval, A, Jd, Js) = Jd / Js                     # the weight scales with the interval

# ---------------------------------------------------------------------------------------

"""
    subdivide(r, domain, n)
    subdivide(r, n)

Composite rule: `domain` split into `n` congruent pieces per direction (`n` intervals,
`n²` triangles), with `r` mapped affinely onto each. The claim is preserved.
"""
subdivide(r::QuadratureRule, n::Integer) = subdivide(r, r.domain, n)

function subdivide(r::QuadratureRule{1,T,<:Interval}, dom::Interval, n::Integer) where {T}
    n >= 1 || throw(ArgumentError("n must be positive"))
    S = _maptype(T, typeof(dom.a))
    return _at_rule_precision(r) do
        h = (S(dom.b) - S(dom.a)) / n
        pieces = [map_to(r, Interval(S(dom.a) + (k - 1) * h, S(dom.a) + k * h)) for k in 1:n]
        _concat(pieces, r, convert_domain(S, dom), "subdivide: $n intervals")
    end
end

function subdivide(r::QuadratureRule{2,T,<:Simplex}, dom::Simplex{2}, n::Integer) where {T}
    n >= 1 || throw(ArgumentError("n must be positive"))
    S = _maptype(T, eltype(eltype(dom.vertices)))
    return _at_rule_precision(r) do
        v0, v1, v2 = (SVector{2,S}(v) for v in dom.vertices)
        P(i, j) = v0 + (S(i) / n) * (v1 - v0) + (S(j) / n) * (v2 - v0)
        cells = Simplex{2,S,3}[]
        for j in 0:(n - 1), i in 0:(n - 1 - j)
            push!(cells, Simplex(P(i, j), P(i + 1, j), P(i, j + 1)))
            i + j <= n - 2 && push!(cells, Simplex(P(i + 1, j + 1), P(i, j + 1), P(i + 1, j)))
        end
        pieces = [map_to(r, c) for c in cells]
        _concat(pieces, r, convert_domain(S, dom), "subdivide: $(n^2) triangles")
    end
end

function subdivide(r::QuadratureRule{3,T,<:Simplex}, dom::Simplex{3}, n::Integer) where {T}
    ispow2(n) || throw(ArgumentError("a tetrahedron is subdivided by repeated bisection, so n must be a " *
                                     "power of two (got $n); n = $(nextpow(2, n)) would give $(nextpow(2, n)^3) cells"))
    S = _maptype(T, eltype(eltype(dom.vertices)))
    return _at_rule_precision(r) do
        cells = [convert_domain(S, dom)]
        for _ in 1:Int(log2(n))
            cells = reduce(vcat, map(bisect_tetrahedron, cells))
        end
        pieces = [map_to(r, c) for c in cells]
        _concat(pieces, r, convert_domain(S, dom), "subdivide: $(n^3) tetrahedra")
    end
end

subdivide(r::QuadratureRule, dom::WeightedDomain, n::Integer) =
    throw(ArgumentError("a rule for a weight cannot be subdivided: the weight is not the same function on each " *
                        "piece, so the pieces need rules of their own; build one per piece with `rule`"))
subdivide(r::QuadratureRule, dom::Domain, n::Integer) =
    throw(ArgumentError("subdivide splits intervals, triangles and tetrahedra; not $(dom)"))

"""
    bisect_tetrahedron(t) -> 8 tetrahedra

Freudenthal bisection: the four corner tetrahedra, plus the inner octahedron split along the
`m₁₂–m₃₄` diagonal. Every child is similar in size to `t/2`, so repeated bisection does not
degenerate.
"""
function bisect_tetrahedron(t::Simplex{3})
    v = t.vertices
    m(i, j) = (v[i] + v[j]) / 2
    m12, m13, m14, m23, m24, m34 = m(1, 2), m(1, 3), m(1, 4), m(2, 3), m(2, 4), m(3, 4)
    return [Simplex(v[1], m12, m13, m14), Simplex(v[2], m12, m23, m24),
            Simplex(v[3], m13, m23, m34), Simplex(v[4], m14, m24, m34),
            Simplex(m12, m13, m14, m34), Simplex(m12, m14, m24, m34),
            Simplex(m12, m24, m23, m34), Simplex(m12, m23, m13, m34)]
end

function _concat(pieces, r, dom, step)
    xs = reduce(vcat, [p.nodes for p in pieces])
    ws = reduce(vcat, [p.weights for p in pieces])
    return QuadratureRule(xs, ws, dom, r.exactness, with_step(r.provenance, step), r.certificate)
end

# ---------------------------------------------------------------------------------------

"""
    transform(r, φ, Jφ; domain = domain(r), claim = nothing)

General change of variables: nodes `φ(x)`, weights `w · |Jφ(x)|`, where `Jφ` returns the
Jacobian matrix or its determinant. The result has **no exactness claim**: under a
nonlinear φ the rule is exact on `{p ∘ φ⁻¹ · |J_{φ⁻¹}|}`, which is not a polynomial space.
Pass `claim` to assert one anyway; the assertion is recorded in the provenance as an
assertion, not a derivation.
"""
function transform(r::QuadratureRule{D}, φ, Jφ; domain::Domain{D} = r.domain, claim = nothing) where {D}
    newclaim = claim === nothing ? NoClaim() : claim
    step = claim === nothing ? "transform: nonlinear change of variables (claim destroyed → NoClaim)" :
           "transform: nonlinear change of variables; claim $(describe(claim)) ASSERTED by the caller, not derived"
    # φ and Jφ run at the rule's precision, so what they compute in BigFloat is not rounded to
    # the ambient one
    return _at_rule_precision(r) do
        xs = [φ(x) for x in r.nodes]
        ws = [w * abs(_detval(Jφ(x))) for (x, w) in zip(r.nodes, r.weights)]
        T = promote_type(eltype(eltype(xs)), eltype(ws))
        xsT = D == 1 ? T.(xs) : [SVector{D,T}(x) for x in xs]
        QuadratureRule(xsT, T.(ws), domain, newclaim, with_step(r.provenance, step; symmetry = NoSymmetry()), r.certificate)
    end
end
_detval(J::Number) = J
_detval(J::AbstractMatrix) = det(J)

"""
    duffy(r; power = 2)

Radial (Duffy-type) grading of a rule on a reference simplex towards the vertex at the
origin: `x ↦ u^{q-1} x` with `u = Σx`, weights `w · q u^{(q-1)D}`. Clusters nodes at the
vertex, which cancels an `r^{-(q-1)D}`-type singularity there. The result has
**no exactness claim**.
"""
function duffy(r::QuadratureRule{D,T,<:Simplex}; power::Integer = 2) where {D,T}
    isreference(r.domain) || throw(ArgumentError("duffy acts on rules on the reference simplex; map afterwards"))
    power >= 1 || throw(ArgumentError("power must be ≥ 1"))
    q = power
    return _at_rule_precision(r) do
        xs = SVector{D,T}[]
        ws = T[]
        for (x, w) in zip(r.nodes, r.weights)
            u = sum(x)
            push!(xs, u^(q - 1) * x)
            push!(ws, w * q * u^((q - 1) * D))
        end
        QuadratureRule(xs, ws, r.domain, NoClaim(),
                       with_step(r.provenance, "duffy: radial grading x ↦ u^$(q-1)·x towards vertex 0 (claim destroyed → NoClaim)";
                                 symmetry = NoSymmetry()), r.certificate)
    end
end
