"""
    CubatureRules

Quadrature and cubature rules generated on demand, at arbitrary order and arbitrary
precision, through one pipeline: seed → refine → certify.

```julia
using CubatureRules
r = rule(Simplex{2}(); degree = 20)
integrate(x -> exp(x[1] * x[2]), r)
```
"""
module CubatureRules

using LinearAlgebra
using Printf
using Random
using SHA: sha256
using StaticArrays
import GenericLinearAlgebra
import InteractiveUtils
import SpecialFunctions
import TOML
import QuadGK
import Test

# core types
include("core/claims.jl")
include("domains/domains.jl")
include("domains/moments.jl")
include("domains/orthobasis.jl")
include("domains/orthobasis3.jl")
include("core/records.jl")
include("core/rule.jl")
include("core/precision.jl")
include("registry/interface.jl")

# refinement
include("refine/newton.jl")

# families and the symmetry machinery they need
include("families/onedim/gaussjacobi.jl")
include("families/onedim/newtoncotes.jl")
include("families/onedim/fejer.jl")
include("families/onedim/delegated.jl")
include("families/simplex/conical.jl")
include("symmetry/orbits.jl")
include("symmetry/invariant.jl")
include("symmetry/moment_system.jl")
include("refine/seeds.jl")
include("refine/elimination.jl")
include("families/simplex/grundmannmoller.jl")
include("families/simplex/xiaogimbutas.jl")
include("families/simplex/fullysymmetric.jl")
include("composition/tensor.jl")

# selection, application, verification, presentation
include("registry/registry.jl")
include("apply/transport.jl")
include("apply/integrate.jl")
include("apply/sequence.jl")
include("verify/verify.jl")
include("emit/show.jl")
include("emit/cite.jl")
include("benchmark/construction.jl")

# domains and claims
export Domain, Interval, Simplex, WeightedDomain, JacobiWeight,
       Orthotope, Sphere, Ball, Polytope, Wedge, Pyramid
export ExactnessClaim, PolynomialDegree, SpanOf, NoClaim
export measure, vertices, barycentric, cartesian, indomain, isinterior, isreference,
       monomial_moment, barycentric_moment, AffineMap, affine_map
# rules
export QuadratureRule, StaticQuadratureRule, static, nodes, weights, domain, exactness,
       provenance, certificate, npoints, degree, family, derivation, rule_hash
export Provenance, Certificate, Verification, Citation, Derived, Seeded
# families and the registry
export RuleFamily, CombinatorFamily, GaussJacobi, GaussLegendre, ConicalProduct, TensorProduct, ⊗,
       NewtonCotes, Fejer, Lobatto, Radau, ClenshawCurtis, GaussKronrod,
       GrundmannMöller, GrundmannMoeller, XiaoGimbutas, FullySymmetric
export rule, available, compare, candidates, properties, degree_range, cost_estimate,
       families
export CancellationToken, cancel!, CancelledError, RefinementError, NoRuleError
export SeedSource, TableSeed, MultistartSeed, ExplicitSeed, LowerDegreeSeed
# application and transport
export integrate, map_to, subdivide, transform, duffy, RuleSequence, EmbeddedRule, embedded,
       IntegrationResult
# verification and presentation
export verify, check, passed, @test_exact, cite
export benchmark_construction

end # module
