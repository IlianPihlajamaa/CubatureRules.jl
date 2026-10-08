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
import InteractiveUtils
import SpecialFunctions
import TOML
import QuadGK
import Test

# core types
include("core/claims.jl")
include("domains/domains.jl")
include("domains/momentweight.jl")
include("domains/functionweight.jl")
include("domains/kernel.jl")
include("domains/surface_triangle.jl")
include("domains/inverse_distance.jl")
include("domains/oscillatory.jl")
include("domains/singular.jl")
include("domains/moments.jl")
include("core/mpfr.jl")
include("domains/orthobasis.jl")
include("domains/orthobasis3.jl")
include("domains/unbounded.jl")
include("domains/sphere.jl")
include("domains/harmonics.jl")
include("domains/ball.jl")
include("domains/wedge_pyramid.jl")
include("core/records.jl")
include("core/rule.jl")
include("core/precision.jl")
include("registry/interface.jl")

# refinement
include("refine/newton.jl")

# families and the symmetry machinery they need
include("families/onedim/gauss_core.jl")
include("families/onedim/gaussjacobi.jl")
include("families/onedim/modified_chebyshev.jl")
include("families/onedim/stieltjes.jl")
include("families/onedim/singular_gauss.jl")
include("families/onedim/filon.jl")
include("families/onedim/newtoncotes.jl")
include("families/onedim/fejer.jl")
include("families/onedim/patterson.jl")
include("families/onedim/tanhsinh.jl")
include("families/onedim/laguerre_hermite.jl")
include("families/onedim/christoffel.jl")
include("families/onedim/desinh.jl")
include("families/onedim/delegated.jl")
include("families/simplex/conical.jl")
include("families/simplex/duffy_gauss.jl")
include("families/simplex/duffy_sinh.jl")
include("symmetry/orbits.jl")
include("symmetry/box.jl")
include("symmetry/invariant.jl")
include("symmetry/moment_system.jl")
include("symmetry/octahedral.jl")
include("symmetry/octahedral_harmonics.jl")
include("symmetry/round.jl")
include("symmetry/pyramid_orbits.jl")
include("symmetry/wedge_orbits.jl")
include("symmetry/octahedral_classify.jl")
include("refine/seeds.jl")
include("refine/elimination.jl")
include("refine/box_search.jl")
include("refine/orbit_search.jl")
include("refine/octahedral_search.jl")
include("refine/octahedral_grow.jl")
include("families/simplex/grundmannmoller.jl")
include("families/simplex/xiaogimbutas.jl")
include("families/simplex/fullysymmetric.jl")
include("families/box/fullysymmetric_box.jl")
include("families/box/fullysymmetric_round.jl")
include("families/wedge_pyramid/fullysymmetric_pyramid.jl")
include("families/wedge_pyramid/fullysymmetric_wedge.jl")
include("families/sphere/product.jl")
include("families/sphere/lebedev.jl")
include("families/sphere/ballproduct.jl")
include("families/sphere/gaussianproduct.jl")
include("families/sphere/exponentialproduct.jl")
include("composition/tensor.jl")
include("composition/wedge.jl")

# selection, application, verification, presentation
include("registry/registry.jl")
include("apply/transport.jl")
include("apply/integrate.jl")
include("apply/embedded.jl")
include("verify/verify.jl")
include("emit/show.jl")
include("emit/cite.jl")
include("benchmark/construction.jl")

# domains and claims
# What `using CubatureRules` brings into scope: the names an ordinary call needs. The rest of
# the public API is declared `public` below, reachable as `CubatureRules.name` or through an
# explicit `using CubatureRules: name`, without crowding the caller's namespace. Quadrature
# and geometry packages collide readily over `vertices`, `⊗` and the like, and a name that is
# rarely typed is not worth a collision.
export Domain, Interval, Simplex, Orthotope, Sphere, Ball, Disk, WeightedDomain, JacobiWeight,
       MomentWeight, OrdinaryMoments, LogWeight, FunctionWeight, PointMass, PrincipalValue, FinitePart, Oscillatory,
       InverseDistance, SurfaceTriangle,
       HalfLine, RealLine, RealSpace, LaguerreRay, HermiteLine, GaussianSpace, ExponentialSpace,
       Polytope, Wedge, Pyramid
export ExactnessClaim, PolynomialDegree, SpanOf, NoClaim
export measure, indomain, isinterior
# rules
export QuadratureRule, nodes, weights, domain, exactness, provenance, certificate,
       npoints, degree, family
# families and the registry
export RuleFamily, GaussJacobi, GaussLegendre, ConicalProduct, TensorProduct,
       NewtonCotes, Fejer, GaussPatterson, TanhSinh, Lobatto, Radau, ClenshawCurtis, GaussKronrod,
       GaussLaguerre, GaussHermite, ExpSinh, SinhSinh,
       GrundmannMöller, GrundmannMoeller, XiaoGimbutas, FullySymmetric, SphereProduct,
       LebedevRule, UpstreamLebedev, BallProduct, GaussianProduct, ExponentialProduct, ModifiedChebyshev,
       StieltjesDiscretization, SingularGauss, Filon, WedgeProduct, DuffyGauss, DuffySinh
export rule, available, compare, families
export RefinementError, NoRuleError
# application and transport
export integrate, map_to, subdivide, transform, duffy,
       EmbeddedRule, embedded, IntegrationResult
# verification and presentation
export verify, check, passed, @test_exact, cite

# Public, but not exported: documented API that most callers never type. `nodes` and
# `weights` stay exported despite colliding with other quadrature packages — they are what
# every caller touches, and that collision is inherent rather than clutter.
public MonicRecurrence, monic, monomial_recurrence, shift, shifted_legendre_recurrence, legendre_monic,
       christoffel, second_kind, endpoints, wheeler, MomentDomain, MomentBreakdownError,
       FunctionDomain, WeightPiece, discretize, discrete_stieltjes, verification_weight,
       SingularKernel, KernelDomain, jacobi_hilbert, kernel_moments, InverseDistanceDomain,
       OscillatoryWeight, OscillatoryDomain, oscillatory_moments, spherical_bessel_up, spherical_bessel_down,
       Recurrence, jacobi_recurrence, laguerre_recurrence, hermite_recurrence,
       ExponentialWeight, GaussianWeight, RadialExponentialWeight, vertices, barycentric, cartesian, isreference,
       monomial_moment, barycentric_moment, AffineMap, affine_map, affine_frame,
       rule_hash,
       Provenance, Certificate, Verification, Citation, Derived, Seeded,
       CombinatorFamily, ⊗,
       candidates, properties, degree_range, cost_estimate, selectable, family_license,
       selection_warnings!, unloaded_alternatives, weight_amplification, POSITIVE_PREFERENCE,
       derivation,
       # the rest of the family interface: what docs/src/tutorial/families.md tells authors to use
       build, BuildContext, outtype, isexact, finalize_number, checkcancel, claimed_degree,
       degree_for_npoints, supports_type, describe_family, reference, verification_basis,
       verify_convergence, classify_octahedral,
       CancellationToken, cancel!, CancelledError,
       SeedSource, TableSeed, MultistartSeed, ExplicitSeed, LowerDegreeSeed,
       benchmark_construction

# A family list cached while precompiling carries a world age from another session.
__init__() = reset_families_cache()

include("precompile.jl")

end # module
