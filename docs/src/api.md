# API

```@meta
CurrentModule = CubatureRules
```

This page documents every exported and public name, grouped by what it is used for. The
rule families and the domain types are documented with their domain in the
[Catalogue](catalogue/index.md).

Names that are `public` but not exported can be used as `CubatureRules.name`, or imported
with `using CubatureRules: name`.

```@docs
CubatureRules
```

## Getting a rule

```@docs
rule
available
compare
families
```

## Rules

```@docs
QuadratureRule
nodes
weights
domain
npoints
degree
exactness
family
provenance
certificate
rule_hash
```

## Exactness claims

```@docs
ExactnessClaim
PolynomialDegree
SpanOf
NoClaim
```

## Integration

```@docs
integrate
IntegrationResult
EmbeddedRule
embedded
```

## Mapping and transformations

```@docs
map_to
subdivide
transform
duffy
⊗
AffineMap
affine_map
```

## Domains

The concrete domains are documented in the [Catalogue](catalogue/index.md).

```@docs
Domain
measure
indomain
isinterior
isreference
reference
vertices
barycentric
cartesian
ExponentialWeight
GaussianWeight
RadialExponentialWeight
MomentDomain
Polytope
Wedge
Pyramid
```

## Measures given by moments

```@docs
MonicRecurrence
monic
shift
shifted_legendre_recurrence
legendre_monic
second_kind
endpoints
monomial_recurrence
wheeler
MomentBreakdownError
Recurrence
jacobi_recurrence
laguerre_recurrence
hermite_recurrence
```

## Verification

```@docs
check
verify
passed
@test_exact
verify_convergence
verification_basis
Verification
monomial_moment
barycentric_moment
```

## Provenance and citation

```@docs
Provenance
Certificate
Citation
cite
license_warnings!
```

## Controlling construction

```@docs
CancellationToken
cancel!
CancelledError
benchmark_construction
```

## Errors

```@docs
NoRuleError
RefinementError
```

## Defining a family

These are the methods a new family implements, and the helpers it uses; see
[I want to build a family](tutorial/families.md).

```@docs
RuleFamily
CombinatorFamily
candidates
build
properties
degree_range
derivation
Derived
Seeded
claimed_degree
degree_for_npoints
cost_estimate
supports_type
describe_family
selectable
family_license
BuildContext
outtype
isexact
finalize_number
checkcancel
classify_octahedral
```

## Seed sources

```@docs
SeedSource
TableSeed
ExplicitSeed
MultistartSeed
LowerDegreeSeed
```
