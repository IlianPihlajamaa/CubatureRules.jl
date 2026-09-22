# Transporting and combining rules, each with a stated effect on the claim (PLAN §2.3):
#
#     map_to(r, dom)        affine only          claim preserved
#     subdivide(r, dom, n)  affine pieces        claim preserved (per cell, hence globally)
#     transform(r, φ, Jφ)   general φ            → NoClaim, unless the caller asserts one
#     duffy(r)              radial grading       → NoClaim

_rule_bits(r::QuadratureRule{D,BigFloat}) where {D} = precision(first(r.weights))
_rule_bits(r::QuadratureRule) = 256

"Run `f` at the precision of the rule's BigFloats (so no conversion reads ambient precision)."
_at_rule_precision(f, r::QuadratureRule) = eltype(r) === BigFloat ? with_bits(f, _rule_bits(r)) : f()

"""
    map_to(r, domain)

The affine image of rule `r` on `domain` (same domain kind and dimension), weights scaled
by `|det A|`. The claim is preserved. There is deliberately no way to pass a nonlinear map
here: use [`transform`](@ref), which reports `NoClaim`.
"""
function map_to(r::QuadratureRule{D,T,<:Simplex}, dom::Simplex{D}) where {D,T}
    S = promote_type(T, eltype(eltype(dom.vertices)))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        A, b, J = _affine(r.domain, dom, S)
        xs = [A * SVector{D,S}(x) + b for x in r.nodes]
        ws = [S(w) * J for w in r.weights]
        QuadratureRule(xs, ws, convert_domain(S, dom), r.exactness,
                       with_step(r.provenance, "map_to: affine image onto $(dom)"), r.certificate)
    end
end

function map_to(r::QuadratureRule{1,T,<:Interval}, dom::Interval) where {T}
    S = promote_type(T, typeof(dom.a))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        s = (S(dom.b) - S(dom.a)) / (S(r.domain.b) - S(r.domain.a))
        xs = [S(dom.a) + s * (S(x) - S(r.domain.a)) for x in r.nodes]
        ws = [S(w) * s for w in r.weights]
        QuadratureRule(xs, ws, convert_domain(S, dom), r.exactness,
                       with_step(r.provenance, "map_to: affine image onto $(dom)"), r.certificate)
    end
end

map_to(r::QuadratureRule, dom::Domain) =
    throw(ArgumentError("map_to transports a rule affinely between domains of the same kind; " *
                        "cannot map a rule on $(r.domain) onto $(dom)"))

# (A, b, |det A|) of the affine map from `src` onto `dst`, in type S. Allocation-free for
# isbits S; the reference source needs no inversion.
@inline function _affine(src::Simplex{D}, dst::Simplex{D}, ::Type{S}) where {D,S}
    Ed = SMatrix{D,D,S}(edge_matrix(dst))
    if isreference(src)
        A = Ed
        b = SVector{D,S}(dst.vertices[1])
    else
        Es = SMatrix{D,D,S}(edge_matrix(src))
        A = Ed / Es
        b = SVector{D,S}(dst.vertices[1]) - A * SVector{D,S}(src.vertices[1])
    end
    return A, b, abs(det(A))
end

@inline function _affine(src::Interval, dst::Interval, ::Type{S}) where {S}
    s = (S(dst.b) - S(dst.a)) / (S(src.b) - S(src.a))
    return s, S(dst.a) - s * S(src.a), abs(s)
end

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
    S = promote_type(T, typeof(dom.a))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        h = (S(dom.b) - S(dom.a)) / n
        pieces = [map_to(r, Interval(S(dom.a) + (k - 1) * h, S(dom.a) + k * h)) for k in 1:n]
        _concat(pieces, r, convert_domain(S, dom), "subdivide: $n intervals")
    end
end

function subdivide(r::QuadratureRule{2,T,<:Simplex}, dom::Simplex{2}, n::Integer) where {T}
    n >= 1 || throw(ArgumentError("n must be positive"))
    S = promote_type(T, eltype(eltype(dom.vertices)))
    S <: Integer && (S = Rational{BigInt})
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
    S = promote_type(T, eltype(eltype(dom.vertices)))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        cells = [convert_domain(S, dom)]
        for _ in 1:Int(log2(n))
            cells = reduce(vcat, map(bisect_tetrahedron, cells))
        end
        pieces = [map_to(r, c) for c in cells]
        _concat(pieces, r, convert_domain(S, dom), "subdivide: $(n^3) tetrahedra")
    end
end

subdivide(r::QuadratureRule, dom::Domain, n::Integer) =
    throw(NotYetImplemented("subdivide for $(typeof(dom).name.name) domains of this dimension", "a later release"))

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
    xs = [D == 1 ? φ(x) : φ(x) for x in r.nodes]
    ws = [w * abs(_detval(Jφ(x))) for (x, w) in zip(r.nodes, r.weights)]
    T = promote_type(eltype(eltype(xs)), eltype(ws))
    newclaim = claim === nothing ? NoClaim() : claim
    step = claim === nothing ? "transform: nonlinear change of variables (claim destroyed → NoClaim)" :
           "transform: nonlinear change of variables; claim $(describe(claim)) ASSERTED by the caller, not derived"
    xsT = D == 1 ? T.(xs) : [SVector{D,T}(x) for x in xs]
    return QuadratureRule(xsT, T.(ws), domain, newclaim, with_step(r.provenance, step; symmetry = :none), r.certificate)
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
                                 symmetry = :none), r.certificate)
    end
end

@inline function _affine(src::Orthotope{D}, dst::Orthotope{D}, ::Type{S}) where {D,S}
    s = SVector{D,S}((dst.hi - dst.lo) ./ (src.hi - src.lo))
    A = SMatrix{D,D,S}(Diagonal(s))
    return A, SVector{D,S}(dst.lo) - A * SVector{D,S}(src.lo), abs(prod(s))
end

# A sphere maps onto another by a similarity — translate and scale — so the Jacobian is
# r^(D-1), the surface one, and the claim survives: a polynomial of degree d pulled back
# through x ↦ c + ρx is again a polynomial of degree d.
function map_to(r::QuadratureRule{D,T,<:Sphere}, dom::Sphere{D}) where {D,T}
    S = promote_type(T, eltype(dom.centre), typeof(dom.radius))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        src = r.domain
        ρ = S(dom.radius) / S(src.radius)
        c = SVector{D,S}(map(S, dom.centre))
        c0 = SVector{D,S}(map(S, src.centre))
        xs = [c + ρ * (SVector{D,S}(map(S, x)) - c0) for x in r.nodes]
        ws = [S(w) * ρ^(D - 1) for w in r.weights]
        QuadratureRule(xs, ws, convert_domain(S, dom), r.exactness,
                       with_step(r.provenance, "map_to: similarity image onto $(dom)"), r.certificate)
    end
end

function map_to(r::QuadratureRule{D,T,<:Orthotope}, dom::Orthotope{D}) where {D,T}
    S = promote_type(T, eltype(dom.lo))
    S <: Integer && (S = Rational{BigInt})
    return _at_rule_precision(r) do
        A, b, J = _affine(r.domain, dom, S)
        xs = [A * SVector{D,S}(x) + b for x in r.nodes]
        ws = [S(w) * J for w in r.weights]
        QuadratureRule(xs, ws, convert_domain(S, dom), r.exactness,
                       with_step(r.provenance, "map_to: affine image onto $(dom)"), r.certificate)
    end
end
