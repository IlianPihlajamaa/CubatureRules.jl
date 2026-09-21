# Exact monomial moments.
#
# Over the reference D-simplex the Dirichlet integral gives, for barycentric exponents β,
#
#     ∫ λ₀^β₀ ⋯ λ_D^β_D dx = β₀! ⋯ β_D! / (|β| + D)!
#
# which also covers Cartesian monomials x^α (take β₀ = 0). Everything is exact
# Rational{BigInt}; conversion to a floating type happens once, at the caller.

"""
    barycentric_moment(β)

Exact `∫ λ^β dx` over the reference simplex of dimension `length(β) - 1`, as a
`Rational{BigInt}`.
"""
function barycentric_moment(β)
    D = length(β) - 1
    D >= 1 || throw(ArgumentError("need at least two barycentric exponents"))
    all(>=(0), β) || throw(ArgumentError("exponents must be non-negative"))
    num = prod(factorial(big(b)) for b in β)
    return num // factorial(big(sum(β) + D))
end

"""
    monomial_moment(dom, α)

Exact `∫ x^α dx` over the reference domain `dom`, as a `Rational{BigInt}`.
"""
function monomial_moment(s::Simplex{D}, α) where {D}
    isreference(s) || throw(ArgumentError("exact moments are provided on the reference simplex; map the rule instead"))
    length(α) == D || throw(ArgumentError("need $D exponents"))
    return barycentric_moment((0, α...))
end

function monomial_moment(d::Interval, α)
    isreference(d) || throw(ArgumentError("exact moments are provided on the reference interval"))
    k = only(α)
    return isodd(k) ? big(0) // 1 : big(2) // (k + 1)
end

# ---------------------------------------------------------------------------------------
# Affine maps. `map_to` is restricted to these at the type level (PLAN §2.3): there is no
# entry point that takes a user function and still reports a polynomial degree.

"""
    AffineMap(A, b)

The map `x ↦ A x + b`. The only kind of map under which a `PolynomialDegree` claim is
preserved.
"""
struct AffineMap{D,T,L}
    A::SMatrix{D,D,T,L}
    b::SVector{D,T}
end

(m::AffineMap)(x) = m.A * x + m.b
jacobian_det(m::AffineMap) = abs(det(m.A))

"""
    affine_map(from, to)

The affine map taking domain `from` onto domain `to` (vertex `i` to vertex `i`).
"""
function affine_map(from::Simplex{D}, to::Simplex{D}) where {D}
    T = promote_type(eltype(eltype(from.vertices)), eltype(eltype(to.vertices)))
    T = T <: Integer ? Rational{BigInt} : T
    Ef = SMatrix{D,D,T}(edge_matrix(from))
    Et = SMatrix{D,D,T}(edge_matrix(to))
    A = Et / Ef
    b = SVector{D,T}(to.vertices[1]) - A * SVector{D,T}(from.vertices[1])
    return AffineMap(A, b)
end

function affine_map(from::Interval, to::Interval)
    T = promote_type(typeof(from.a), typeof(to.a))
    T = T <: Integer ? Rational{BigInt} : T
    s = T(to.b - to.a) / T(from.b - from.a)
    return AffineMap(SMatrix{1,1,T}(s), SVector{1,T}(T(to.a) - s * T(from.a)))
end

function monomial_moment(d::Orthotope{D}, α) where {D}
    isreference(d) || throw(ArgumentError("exact moments are provided on the reference orthotope"))
    length(α) == D || throw(ArgumentError("need $D exponents"))
    return prod(isodd(k) ? big(0) // 1 : big(2) // (k + 1) for k in α)
end

"The affine map between two boxes is diagonal."
function affine_map(from::Orthotope{D}, to::Orthotope{D}) where {D}
    T = promote_type(eltype(from.lo), eltype(to.lo))
    T = T <: Integer ? Rational{BigInt} : T
    s = SVector{D,T}((to.hi - to.lo) ./ (from.hi - from.lo))
    A = SMatrix{D,D,T}(Diagonal(s))
    return AffineMap(A, SVector{D,T}(to.lo) - A * SVector{D,T}(from.lo))
end
