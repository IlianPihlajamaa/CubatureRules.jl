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
include("domains/unbounded.jl")
include("domains/sphere.jl")
include("domains/harmonics.jl")
include("domains/ball.jl")
include("core/records.jl")
include("core/rule.jl")
include("core/precision.jl")
include("registry/interface.jl")

# refinement
include("refine/newton.jl")

# families and the symmetry machinery they need
include("families/onedim/gauss_core.jl")
include("families/onedim/gaussjacobi.jl")
include("families/onedim/newtoncotes.jl")
include("families/onedim/fejer.jl")
include("families/onedim/tanhsinh.jl")
include("families/onedim/laguerre_hermite.jl")
include("families/onedim/desinh.jl")
include("families/onedim/delegated.jl")
include("families/simplex/conical.jl")
include("symmetry/orbits.jl")
include("symmetry/invariant.jl")
include("symmetry/moment_system.jl")
include("symmetry/octahedral.jl")
include("refine/seeds.jl")
include("refine/elimination.jl")
include("refine/octahedral_search.jl")
include("families/simplex/grundmannmoller.jl")
include("families/simplex/xiaogimbutas.jl")
include("families/simplex/fullysymmetric.jl")
include("families/sphere/product.jl")
include("families/sphere/lebedev.jl")
include("families/sphere/ballproduct.jl")
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
export Domain, Interval, Simplex, WeightedDomain, JacobiWeight, HalfLine, RealLine,
       ExponentialWeight, GaussianWeight, LaguerreRay, HermiteLine,
       Orthotope, Sphere, Ball, Disk, Polytope, Wedge, Pyramid
export ExactnessClaim, PolynomialDegree, SpanOf, NoClaim
export measure, vertices, barycentric, cartesian, indomain, isinterior, isreference,
       monomial_moment, barycentric_moment, AffineMap, affine_map
# rules
export QuadratureRule, StaticQuadratureRule, static, nodes, weights, domain, exactness,
       provenance, certificate, npoints, degree, family, derivation, rule_hash
export Provenance, Certificate, Verification, Citation, Derived, Seeded
# families and the registry
export RuleFamily, CombinatorFamily, GaussJacobi, GaussLegendre, ConicalProduct, TensorProduct, ⊗,
       NewtonCotes, Fejer, TanhSinh, Lobatto, Radau, ClenshawCurtis, GaussKronrod,
       GaussLaguerre, GaussHermite, ExpSinh, SinhSinh,
       GrundmannMöller, GrundmannMoeller, XiaoGimbutas, FullySymmetric, SphereProduct, Lebedev, BallProduct
export rule, available, compare, candidates, properties, degree_range, cost_estimate,
       families
export CancellationToken, cancel!, CancelledError, RefinementError, NoRuleError
export SeedSource, TableSeed, MultistartSeed, ExplicitSeed, LowerDegreeSeed
# application and transport
export integrate, map_to, subdivide, transform, duffy, RuleSequence, LevelSequence, EmbeddedRule, embedded,
       IntegrationResult
# verification and presentation
export verify, check, passed, @test_exact, cite
export benchmark_construction

end # module
