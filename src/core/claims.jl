# What a rule claims (PLAN §2.2).
#
# The claim hierarchy models *exactness* only. Convergence behaviour is a property of a
# family of rules and belongs to a (future) error subsystem, not here.

"""
    ExactnessClaim

What a quadrature rule claims to integrate exactly. Subtypes:

- [`PolynomialDegree`](@ref) — exact on all polynomials of total degree `≤ d` (against the
  domain's measure, which may carry a weight).
- [`SpanOf`](@ref) — exact on the span of a given finite basis.
- [`NoClaim`](@ref) — the rule converges, but is exact on nothing in particular.
"""
abstract type ExactnessClaim end

"""
    PolynomialDegree(d)

Exact on every polynomial of total degree `≤ d` against the measure of the rule's domain.
"""
struct PolynomialDegree <: ExactnessClaim
    d::Int
    function PolynomialDegree(d::Integer)
        d >= 0 || throw(ArgumentError("a polynomial degree must be non-negative, got $d"))
        return new(Int(d))
    end
end

"""
    SpanOf(basis)

Exact on the span of `basis`, a finite collection of functions. The only claim that is
not reducible to a polynomial degree against some measure (e.g. generalised Gaussian
quadrature on Chebyshev systems).
"""
struct SpanOf{B} <: ExactnessClaim
    basis::B
end

"""
    NoClaim()

No exactness claim. Produced by nonlinear transports ([`transform`](@ref), [`duffy`](@ref))
and by families characterised by a convergence rate rather than exactness.
"""
struct NoClaim <: ExactnessClaim end

Base.:(==)(a::SpanOf, b::SpanOf) = a.basis == b.basis
Base.hash(c::SpanOf, h::UInt) = hash(c.basis, hash(:SpanOf, h))

"""
    degree(rule)
    degree(claim)

The polynomial degree of a rule whose claim is a [`PolynomialDegree`](@ref). Errors for
any other claim, since there is no degree to report.
"""
degree(c::PolynomialDegree) = c.d
function degree(c::ExactnessClaim)
    throw(ArgumentError("the claim $(c) has no polynomial degree; `degree` is defined only for " *
                        "PolynomialDegree claims. Inspect `exactness(rule)` instead."))
end

describe(c::PolynomialDegree) = "polynomial degree $(c.d)"
describe(c::SpanOf) = "span of a $(length(c.basis))-function basis"
describe(::NoClaim) = "no exactness claim"
