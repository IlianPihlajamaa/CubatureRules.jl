# A downstream package that adds a rule family. It is precompiled separately from
# CubatureRules, and the test suite asserts that the selector finds it with no
# registration step (PLAN §2.5).
module TestFamilyPackage

using CubatureRules
using StaticArrays
# The family interface is public rather than exported, so it is imported by name — which is
# also what a real downstream package does, and what `docs/src/extending.md` recommends.
import CubatureRules: candidates, build, npoints, properties, degree_range, derivation,
                      isreference, Derived, Provenance

"The one-point centroid rule on any reference simplex: degree 1."
struct CentroidRule <: RuleFamily end

candidates(::Type{CentroidRule}, dom::Simplex, c::PolynomialDegree) =
    isreference(dom) && c.d <= 1 ? [CentroidRule()] : CentroidRule[]
npoints(::CentroidRule, dom, degree) = 1
properties(::CentroidRule, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)
degree_range(::CentroidRule, dom) = 0:1
derivation(::Type{CentroidRule}) = Derived()

function build(::CentroidRule, dom::Simplex{D}, degree::Int, ctx) where {D}
    T = CubatureRules.outtype(ctx)
    x = SVector{D,T}(ntuple(_ -> one(T) / (D + 1), D))
    w = one(T) / factorial(D)
    prov = Provenance(family = "CentroidRule", derivation = Derived())
    return QuadratureRule([x], [w], Simplex{D}(), PolynomialDegree(1), prov)
end

end
