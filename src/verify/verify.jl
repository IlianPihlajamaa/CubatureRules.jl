# Verification (PLAN §8): is this rule what it claims to be? Distinct from the Certificate,
# which only says that Newton converged on the system it was given.
#
# Exactness is tested against *orthonormal* bases (Dubiner on triangles and tetrahedra,
# Legendre/Jacobi on intervals) at twice the rule's precision; monomials would produce false failures at
# high degree. Each basis function's residual is held to a tolerance derived from the
# rounding of the delivered nodes and weights: 16 ε Σ |wᵢ| (|φ(xᵢ)| + |∇φ(xᵢ)|·max(|xᵢ|, 1)).
# Exact
# (Rational) rules are held to zero.

"""
    verify(rule; degree, bits, integrals) -> Verification

Check `rule` against its exactness claim, dispatching on the claim type:

- `PolynomialDegree(d)`: exactness on an orthonormal basis to degree `d`, sharpness (not
  exact at `d + 1`), weights summing to the measure, interior nodes, positive weights, and
  the claimed symmetry group preserved.
- `SpanOf(basis)`: exactness on the given basis; pass `integrals` (the exact integrals).
- `NoClaim`: nothing to verify by exact integration — see [`verify_convergence`](@ref).

`degree` overrides the degree tested (sharpness is then not tested unless it equals the
claim); `bits` overrides the arithmetic precision (default twice the rule's).

A rule on `Sphere{3}` that claims octahedral symmetry is checked on its orbits: the nodes are
grouped into `O_h` orbits (which checks the symmetry), and exactness is tested with
group-averaged even monomials evaluated once per orbit. This is equivalent to the general
check for a symmetric rule and much faster for large ones. `use_symmetry = false` forces the
general check against all spherical harmonics at every node.
"""
verify(r::QuadratureRule; kw...) = verify(r, r.exactness; kw...)

"""
    check(rule; degree, bits, use_symmetry = true) -> Verification

Re-run verification on demand, optionally at a user-specified degree and precision. The
public entry point for rules constructed by hand or loaded from elsewhere. See
[`verify`](@ref) for `use_symmetry`.
"""
check(r::QuadratureRule; degree = nothing, bits = nothing, use_symmetry::Bool = true) =
    verify(r; degree, bits, use_symmetry)

verify(r::QuadratureRule, ::NoClaim; kw...) =
    throw(ArgumentError("a rule with NoClaim has nothing to verify by exact integration; " *
                        "use verify_convergence on a sequence of such rules"))

_is_exact_type(::Type{T}) where {T} = T <: Rational || T <: Integer
_rule_precision(r::QuadratureRule) = _is_exact_type(eltype(r)) ? 0 :
                                     eltype(r) === BigFloat ? precision(first(r.weights)) : bits_of(eltype(r))

# Pointwise basis evaluation on reference domains: returns (φ, |∇φ|₁) vectors for degree
# ≤ n, the exact integrals of the basis, and the index ranges of each degree block.
abstract type VerificationBasis end

struct DubinerBasis{S} <: VerificationBasis
    ws::DubinerWorkspace{S}
    φ::Vector{S}
    gx::Vector{S}
    gy::Vector{S}
end
DubinerBasis{S}(n; normalize) where {S} =
    DubinerBasis{S}(DubinerWorkspace{S}(n; normalize), zeros(S, dubiner_length(n)), zeros(S, dubiner_length(n)),
                    zeros(S, dubiner_length(n)))
function evaluate!(b::DubinerBasis, x)
    dubiner!(b.φ, b.gx, b.gy, b.ws, x[1], x[2])
    return b.φ, abs.(b.gx) .+ abs.(b.gy)
end
blocks(b::DubinerBasis) = [dubiner_index(0, k):dubiner_index(k, 0) for k in 0:(b.ws.n)]

struct JacobiBasis{S,P} <: VerificationBasis
    n::Int
    α::P
    β::P
    normalize::Bool
end
function evaluate!(b::JacobiBasis{S}, xv) where {S}
    x = xv[1]
    p = zeros(S, b.n + 1)
    dp = zeros(S, b.n + 1)
    if b.normalize
        μ = S(jacobi_mass(b.α, b.β))
        A, B = S(b.α), S(b.β)
        p[1] = one(S) / sqrt(μ)
        if b.n >= 1
            b1 = jacobi_b(1, A, B)
            p[2] = (x - jacobi_a(0, A, B)) * p[1] / b1
            dp[2] = p[1] / b1
        end
        for k in 1:(b.n - 1)
            bk, bk1, ak = jacobi_b(k, A, B), jacobi_b(k + 1, A, B), jacobi_a(k, A, B)
            p[k + 2] = ((x - ak) * p[k + 1] - bk * p[k]) / bk1
            dp[k + 2] = (p[k + 1] + (x - ak) * dp[k + 1] - bk * dp[k]) / bk1
        end
    else   # exact arithmetic: unnormalised Legendre, α = β = 0
        p[1] = one(S)
        b.n >= 1 && (p[2] = S(x); dp[2] = one(S))
        for k in 1:(b.n - 1)
            p[k + 2] = ((2k + 1) * x * p[k + 1] - k * p[k]) / (k + 1)
            dp[k + 2] = ((2k + 1) * (p[k + 1] + x * dp[k + 1]) - k * dp[k]) / (k + 1)
        end
    end
    return p, abs.(dp)
end
blocks(b::JacobiBasis) = [k:k for k in 1:(b.n + 1)]

struct BarycentricMonomialBasis{S} <: VerificationBasis
    D::Int
    degrees::Vector{Int}       # homogeneous degrees tested (d, and d+1 for sharpness)
    exps::Vector{Vector{Int}}
    blockranges::Vector{UnitRange{Int}}
end
function BarycentricMonomialBasis{S}(D, degrees) where {S}
    exps = Vector{Vector{Int}}()
    ranges = UnitRange{Int}[]
    for k in degrees
        lo = length(exps) + 1
        append!(exps, compositions(k, D + 1))
        push!(ranges, lo:length(exps))
    end
    return BarycentricMonomialBasis{S}(D, collect(degrees), exps, ranges)
end
function evaluate!(b::BarycentricMonomialBasis{S}, x) where {S}
    λ = vcat(one(S) - sum(x), x)
    φ = [prod(λ[i]^β[i] for i in eachindex(β)) for β in b.exps]
    g = [sum((β[i] * prod(λ[j]^(j == i ? β[j] - 1 : β[j]) for j in eachindex(β)) for i in eachindex(β) if β[i] > 0);
             init = zero(S)) for β in b.exps]
    return φ, abs.(g)
end
blocks(b::BarycentricMonomialBasis) = b.blockranges

struct TetVerifyBasis{S,B<:SimplexBasis{3,S}} <: VerificationBasis
    basis::B
end
TetVerifyBasis{S}(n; normalize) where {S} = (b = SimplexBasis{3,S}(n; normalize); TetVerifyBasis{S,typeof(b)}(b))
function evaluate!(b::TetVerifyBasis, x)
    φ, G = evaluate!(b.basis, x)
    return φ, [abs(G[k, 1]) + abs(G[k, 2]) + abs(G[k, 3]) for k in axes(G, 1)]
end
blocks(b::TetVerifyBasis) = degree_blocks(b.basis)
exact_integrals(b::TetVerifyBasis{S}) where {S} =
    (v = zeros(S, basis_length(b.basis)); v[1] = b.basis.ws.c[1] * S(1 // 6); v)

# Tensor Legendre on the reference box, graded by *total* degree — the space a tensor rule
# claims (PLAN §2.2), not the larger tensor-product space it is exact on.
struct BoxBasis{D,S} <: VerificationBasis
    n::Int
    idx::Vector{NTuple{D,Int}}
    blockranges::Vector{UnitRange{Int}}
    axis::JacobiBasis{S,Int}
end

function BoxBasis{D,S}(n::Int; normalize::Bool = true) where {D,S}
    idx = NTuple{D,Int}[]
    ranges = UnitRange{Int}[]
    for k in 0:n
        lo = length(idx) + 1
        for c in compositions(k, D)
            push!(idx, NTuple{D,Int}(c))
        end
        push!(ranges, lo:length(idx))
    end
    return BoxBasis{D,S}(n, idx, ranges, JacobiBasis{S,Int}(n, 0, 0, normalize))
end

function evaluate!(b::BoxBasis{D,S}, x) where {D,S}
    ps = [evaluate!(b.axis, (x[j],)) for j in 1:D]        # (values, |derivatives|) per axis
    φ = [prod(ps[j][1][i[j] + 1] for j in 1:D) for i in b.idx]
    g = [sum(ps[j][2][i[j] + 1] * prod(abs(ps[l][1][i[l] + 1]) for l in 1:D if l != j; init = one(S))
             for j in 1:D) for i in b.idx]
    return φ, g
end
blocks(b::BoxBasis) = b.blockranges
function exact_integrals(b::BoxBasis{D,S}) where {D,S}
    v = zeros(S, length(b.idx))
    v[1] = (b.axis.normalize ? sqrt(S(2)) : S(2))^D        # only the constant term survives
    return v
end

exact_integrals(b::DubinerBasis{S}) where {S} =
    (v = zeros(S, length(b.φ)); v[1] = b.ws.c[1] * S(1 // 2); v)
function exact_integrals(b::JacobiBasis{S}) where {S}
    v = zeros(S, b.n + 1)
    v[1] = b.normalize ? sqrt(S(jacobi_mass(b.α, b.β))) : S(2)
    return v
end
exact_integrals(b::BarycentricMonomialBasis{S}) where {S} = [S(barycentric_moment(β)) for β in b.exps]

# Orthonormal polynomials of any weight given its recurrence — the basis for verifying
# rules on weighted domains whose polynomials are not Jacobi (Laguerre, Hermite).
struct RecurrenceBasis{S,R<:Recurrence} <: VerificationBasis
    n::Int
    rec::R
end
function evaluate!(b::RecurrenceBasis{S}, xv) where {S}
    x = xv[1]
    p = zeros(S, b.n + 1)
    dp = zeros(S, b.n + 1)
    μ = S(b.rec.mass(S))
    p[1] = one(S) / sqrt(μ)
    if b.n >= 1
        b1 = b.rec.b(1, S)
        p[2] = (x - b.rec.a(0, S)) * p[1] / b1
        dp[2] = p[1] / b1
    end
    for k in 1:(b.n - 1)
        bk, bk1, ak = b.rec.b(k, S), b.rec.b(k + 1, S), b.rec.a(k, S)
        p[k + 2] = ((x - ak) * p[k + 1] - bk * p[k]) / bk1
        dp[k + 2] = (p[k + 1] + (x - ak) * dp[k + 1] - bk * dp[k]) / bk1
    end
    return p, abs.(dp)
end
blocks(b::RecurrenceBasis) = [k:k for k in 1:(b.n + 1)]
exact_integrals(b::RecurrenceBasis{S}) where {S} =
    (v = zeros(S, b.n + 1); v[1] = sqrt(S(b.rec.mass(S))); v)

# Sharpness on an orthonormal basis in one dimension. The test asks whether the rule is exact
# on one function of degree d + 1; every such function has the same error up to its leading
# coefficient, so which one is used matters only for how the error compares with rounding.
# The top basis function p_{d+1} is a poor choice: for an n-point Gauss rule its error is the
# ratio of leading coefficients k_{2n} / k_n², which for Laguerre is (n!)² / (2n)! ≈ 4⁻ⁿ and
# falls below the rounding tolerance once n exceeds about 1.7 times the digits — at which point
# a correct rule was reported as exact one degree too high, and failed. The product
# p_a p_b with a + b = d + 1 has integral 1 or 0 and, for a Gauss rule, is p_n², whose error
# is exactly 1 whatever the weight or n.
struct ProductTop{S,B<:VerificationBasis} <: VerificationBasis
    inner::B
    m::Int                              # the degree tested for sharpness, d + 1
end
ProductTop{S}(inner::B, m::Integer) where {S,B} = ProductTop{S,B}(inner, Int(m))
function evaluate!(b::ProductTop{S}, x) where {S}
    p, dp = evaluate!(b.inner, x)
    φ, g = copy(p), copy(dp)
    i, j = b.m ÷ 2 + 1, b.m - b.m ÷ 2 + 1           # 1-based indices of p_a and p_b
    φ[end] = p[i] * p[j]
    g[end] = dp[i] * abs(p[j]) + abs(p[i]) * dp[j]
    return φ, g
end
blocks(b::ProductTop) = blocks(b.inner)
function exact_integrals(b::ProductTop{S}) where {S}
    v = copy(exact_integrals(b.inner))
    v[end] = iseven(b.m) ? one(S) : zero(S)          # ∫ p_a p_b = δ_ab
    return v
end

# Map a rule to its reference domain in type S: (nodes as Vector{Vector{S}}, weights).
function reference_nodes(r::QuadratureRule{D,T,<:Simplex}, ::Type{S}) where {D,T,S}
    dom = r.domain
    Ed = SMatrix{D,D,S}(edge_matrix(dom))
    v0 = SVector{D,S}(dom.vertices[1])
    J = abs(det(Ed))
    xs = [Vector{S}(Ed \ (SVector{D,S}(map(S, x)) - v0)) for x in r.nodes]
    ws = [S(w) / J for w in r.weights]
    return xs, ws
end
function reference_nodes(r::QuadratureRule{1,T,<:Interval}, ::Type{S}) where {T,S}
    a, b = S(r.domain.a), S(r.domain.b)
    xs = [[(2S(x) - (a + b)) / (b - a)] for x in r.nodes]
    ws = [2S(w) / (b - a) for w in r.weights]
    return xs, ws
end
function reference_nodes(r::QuadratureRule{D,T,<:Orthotope}, ::Type{S}) where {D,T,S}
    lo = SVector{D,S}(map(S, r.domain.lo))
    hi = SVector{D,S}(map(S, r.domain.hi))
    scale = prod((hi - lo) ./ 2)
    xs = [Vector{S}((2 .* SVector{D,S}(map(S, x)) .- (lo .+ hi)) ./ (hi .- lo)) for x in r.nodes]
    return xs, [S(w) / scale for w in r.weights]
end

# A ball maps to its reference the same way, with the volume Jacobian r^D.
function reference_nodes(r::QuadratureRule{D,T,<:Ball}, ::Type{S}) where {D,T,S}
    dom = r.domain
    c = SVector{D,S}(map(S, dom.centre))
    ρ = S(dom.radius)
    xs = [Vector{S}((SVector{D,S}(map(S, x)) .- c) ./ ρ) for x in r.nodes]
    return xs, [S(w) / ρ^D for w in r.weights]
end

# A sphere maps to its reference by translating to the origin and scaling: the weights
# carry r^(D-1), the surface Jacobian, not r^D.
function reference_nodes(r::QuadratureRule{D,T,<:Sphere}, ::Type{S}) where {D,T,S}
    dom = r.domain
    c = SVector{D,S}(map(S, dom.centre))
    ρ = S(dom.radius)
    xs = [Vector{S}((SVector{D,S}(map(S, x)) .- c) ./ ρ) for x in r.nodes]
    return xs, [S(w) / ρ^(D - 1) for w in r.weights]
end

# A weighted domain in more than one dimension is its own reference here: there is no map
# to apply, so the nodes pass through.
function reference_nodes(r::QuadratureRule{D,T,<:WeightedDomain}, ::Type{S}) where {D,T,S}
    isreference(r.domain) || throw(ArgumentError("verification of weighted rules needs the reference domain"))
    return [Vector{S}(map(S, x)) for x in r.nodes], [S(w) for w in r.weights]
end

function reference_nodes(r::QuadratureRule{1,T,<:WeightedDomain}, ::Type{S}) where {T,S}
    isreference(r.domain.base) || throw(ArgumentError("verification of weighted rules needs the reference interval"))
    return [[S(x)] for x in r.nodes], [S(w) for w in r.weights]
end

# A measure given by its moments is stated on the interval it was integrated over, whatever
# that interval is, so the rule is already on its own reference and the nodes pass through.
reference_nodes(r::QuadratureRule{1,T,<:MomentDomain}, ::Type{S}) where {T,S} =
    ([[S(x)] for x in r.nodes], [S(w) for w in r.weights])

"""
    verification_basis(domain, degrees, S, exact) -> (basis, description)

The basis [`verify`](@ref) tests a rule on `domain` against, for the given degrees, in
arithmetic `S`; `exact` is set for rules in exact rational arithmetic. Define a method for a
new kind of domain to make its degree claims verifiable. The description is shown in the
verification report.
"""
function verification_basis end

verification_basis(dom::Simplex{2}, n, ::Type{S}, exact) where {S} = DubinerBasis{S}(n; normalize = !exact),
    exact ? "Dubiner (unnormalised, exact arithmetic)" : "orthonormal Dubiner"
verification_basis(dom::Simplex{3}, n, ::Type{S}, exact) where {S} = TetVerifyBasis{S}(n; normalize = !exact),
    exact ? "tetrahedral Dubiner (unnormalised, exact arithmetic)" : "orthonormal tetrahedral Dubiner"
verification_basis(dom::Simplex{D}, n, ::Type{S}, exact) where {D,S} =
    BarycentricMonomialBasis{S}(D, n), "barycentric monomials of degree $(join(n, ", ")) (exact Dirichlet moments)"
verification_basis(dom::Orthotope{D}, n, ::Type{S}, exact) where {D,S} =
    BoxBasis{D,S}(n; normalize = !exact),
    exact ? "tensor Legendre (unnormalised, exact arithmetic), graded by total degree" :
    "tensor orthonormal Legendre, graded by total degree"
verification_basis(dom::LaguerreDomain, n, ::Type{S}, exact) where {S} =
    RecurrenceBasis{S,typeof(laguerre_recurrence(dom.weight.α))}(n, laguerre_recurrence(dom.weight.α)),
    "orthonormal Laguerre($(dom.weight.α))"
verification_basis(dom::HermiteDomain, n, ::Type{S}, exact) where {S} =
    RecurrenceBasis{S,typeof(hermite_recurrence())}(n, hermite_recurrence()), "orthonormal Hermite"
# Spheres: the harmonics of degree ≤ d span the polynomials of degree ≤ d restricted to the
# sphere, which is what a claim of degree d on a sphere means.
struct HarmonicBasis{S} <: VerificationBasis
    n::Int
    Y::Vector{S}
    P::Vector{S}
    g::Vector{S}                      # the per-degree gradient bound, laid out like Y
end
function HarmonicBasis{S}(n::Integer) where {S}
    g = zeros(S, harmonic_length(n))
    for l in 0:n, i in harmonic_block(l)
        g[i] = harmonic_gradient_bound(l, S)
    end
    return HarmonicBasis{S}(n, zeros(S, harmonic_length(n)), zeros(S, legendre_length(n)), g)
end
function evaluate!(b::HarmonicBasis, x)
    real_harmonics!(b.Y, b.P, b.n, x)
    return b.Y, b.g
end
blocks(b::HarmonicBasis) = [harmonic_block(l) for l in 0:(b.n)]
# Only the constant harmonic integrates to anything: Y₀₀ = 1/√(4π) over an area of 4π.
exact_integrals(b::HarmonicBasis{S}) where {S} =
    (v = zeros(S, harmonic_length(b.n)); v[1] = sqrt(4 * S(π)); v)

struct CircleBasis{S} <: VerificationBasis
    n::Int
    F::Vector{S}
    g::Vector{S}
end
function CircleBasis{S}(n::Integer) where {S}
    g = zeros(S, circle_length(n))
    for k in 0:n, i in circle_block(k)
        g[i] = S(k) / sqrt(S(π))          # |d/dφ| of cos(kφ)/√π, exactly
    end
    return CircleBasis{S}(n, zeros(S, circle_length(n)), g)
end
function evaluate!(b::CircleBasis, x)
    fourier_circle!(b.F, b.n, x)
    return b.F, b.g
end
blocks(b::CircleBasis) = [circle_block(k) for k in 0:(b.n)]
exact_integrals(b::CircleBasis{S}) where {S} =
    (v = zeros(S, circle_length(b.n)); v[1] = sqrt(2 * S(π)); v)

# Cartesian monomials with known moments: the test set for every domain whose orthonormal
# basis either does not exist here or is not worth building. One type serves three domains,
# differing only in which moment function it carries.
#
# On a ball or a Gaussian space these are a genuine basis — the coordinates satisfy no
# relation. On a sphere they are not, since Σxᵢ² = 1 makes them dependent; but a degree claim
# says every polynomial of that degree integrates correctly, and checking that needs a
# spanning test set with known integrals, not an independent one. The redundancy costs
# repeated equations, not correctness.
struct MonomialTestSet{D,S,F} <: VerificationBasis
    exps::Vector{NTuple{D,Int}}
    blockranges::Vector{UnitRange{Int}}
    moment::F
end
function MonomialTestSet{D,S}(degrees, moment::F) where {D,S,F}
    exps = NTuple{D,Int}[]
    ranges = UnitRange{Int}[]
    for k in degrees
        lo = length(exps) + 1
        for c in compositions(k, D)
            push!(exps, NTuple{D,Int}(c))
        end
        push!(ranges, lo:length(exps))
    end
    return MonomialTestSet{D,S,F}(exps, ranges, moment)
end
function evaluate!(b::MonomialTestSet{D,S}, x) where {D,S}
    φ = [prod(S(x[i])^α[i] for i in 1:D) for α in b.exps]
    g = [sum((α[i] * prod(S(x[j])^(j == i ? α[j] - 1 : α[j]) for j in 1:D) for i in 1:D if α[i] > 0);
             init = zero(S)) for α in b.exps]
    return φ, abs.(g)
end
blocks(b::MonomialTestSet) = b.blockranges
exact_integrals(b::MonomialTestSet{D,S}) where {D,S} = [S(b.moment(α)) for α in b.exps]

_degrees(n) = n isa Integer ? (0:n) : n

verification_basis(dom::Ball{D}, n, ::Type{S}, exact) where {D,S} =
    MonomialTestSet{D,S}(_degrees(n), α -> ball_moment(D, α)),
    "monomials of degree $(join(n, ", ")) (exact ball moments)"

verification_basis(dom::GaussianDomain{D}, n, ::Type{S}, exact) where {D,S} =
    MonomialTestSet{D,S}(_degrees(n), α -> gaussian_moment(D, α)),
    "monomials of degree $(join(n, ", ")) (exact Gaussian moments)"

verification_basis(dom::Sphere{3}, n, ::Type{S}, exact) where {S} =
    HarmonicBasis{S}(n), "real spherical harmonics"
verification_basis(dom::Sphere{2}, n, ::Type{S}, exact) where {S} =
    CircleBasis{S}(n), "Fourier basis on the circle"
# Above the 2-sphere there are no hyperspherical harmonics here, and none are needed.
# Monomials cannot be a *basis* on a sphere, since the coordinates satisfy Σxᵢ² = 1, but a
# degree claim says every polynomial of that degree integrates correctly, and checking that
# needs a spanning *test set* with known integrals — which monomials are. The redundancy
# costs repeated equations, not correctness, and the exact moments are already in hand for
# every dimension.
verification_basis(dom::Sphere{D}, n, ::Type{S}, exact) where {D,S} =
    MonomialTestSet{D,S}(_degrees(n), α -> sphere_moment(D, α)),
    "monomials of degree $(join(n, ", ")) (exact sphere moments; dependent, so a test set rather than a basis)"

verification_basis(dom::Interval, n, ::Type{S}, exact) where {S} =
    JacobiBasis{S,Int}(n, 0, 0, !exact), exact ? "Legendre (unnormalised, exact arithmetic)" : "orthonormal Legendre"
verification_basis(dom::WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight}, n, ::Type{S}, exact) where {S} =
    JacobiBasis{S,typeof(dom.weight.α)}(n, dom.weight.α, dom.weight.β, true),
    "orthonormal Jacobi($(dom.weight.α), $(dom.weight.β))"

# A measure known only by its moments is verified against those moments. The auxiliary
# polynomials π₀ … π_d span P_d, so this is the full exactness claim and not a weaker sample
# of it — and the exact integrals are the input data rather than anything derived from it.
struct MomentBasis{S,W<:MomentWeight} <: VerificationBasis
    degs::Vector{Int}
    w::W
end
MomentBasis{S}(degs, w::W) where {S,W<:MomentWeight} = MomentBasis{S,W}(collect(degs), w)
function evaluate!(b::MomentBasis{S}, xv) where {S}
    x = S(xv[1])
    N = maximum(b.degs) + 1
    p = zeros(S, N)
    dp = zeros(S, N)
    p[1] = one(S)
    for k in 1:(N - 1)
        a, c = S(b.w.aux.a(k - 1, S)), S(b.w.aux.b(k - 1, S))
        prev, dprev = k == 1 ? (zero(S), zero(S)) : (p[k - 1], dp[k - 1])
        p[k + 1] = (x - a) * p[k] - c * prev
        dp[k + 1] = p[k] + (x - a) * dp[k] - c * dprev
    end
    return S[p[k + 1] for k in b.degs], S[abs(dp[k + 1]) for k in b.degs]
end
blocks(b::MomentBasis) = [k:k for k in eachindex(b.degs)]
exact_integrals(b::MomentBasis{S}) where {S} = S[S(b.w.moments(k, S)) for k in b.degs]

verification_basis(dom::MomentDomain, n, ::Type{S}, exact) where {S} =
    MomentBasis{S}(_degrees(n), dom.weight),
    "monic auxiliary polynomials of $(dom.weight.label), against the given moments"

"Per-degree-block (max residual, max residual/tolerance, max tolerance)."
function block_residuals(basis, xs, ws, ε, floor_tol)
    I = exact_integrals(basis)
    L = length(I)
    Σ = zeros(eltype(ws), L)
    scale = zeros(eltype(ws), L)
    for (x, w) in zip(xs, ws)
        φ, g = evaluate!(basis, x)
        # a node near zero still carries an absolute rounding error of order ε (it is a
        # difference of O(1) quantities), so the perturbation bound uses max(|x|, 1)
        xnorm = max(maximum(abs, x), one(eltype(x)))
        for k in 1:L
            Σ[k] += w * φ[k]
            scale[k] += abs(w) * (abs(φ[k]) + g[k] * xnorm)
        end
    end
    out = Tuple{BigFloat,BigFloat,BigFloat}[]
    for blk in blocks(basis)
        mres = big(0.0); mratio = big(0.0); mtol = big(0.0)
        for k in blk
            res = BigFloat(abs(Σ[k] - I[k]))
            tol = BigFloat(16ε * scale[k]) + floor_tol
            mres = max(mres, res)
            mtol = max(mtol, tol)
            ratio = iszero(tol) ? (iszero(res) ? big(0.0) : big(Inf)) : res / tol
            mratio = max(mratio, ratio)
        end
        push!(out, (mres, mratio, mtol))
    end
    return out
end

# --- O_h-symmetric rules on the sphere: verification on orbit representatives ---------------
#
# A Lebedev rule of degree 125 has 5294 nodes. Checked the general way, every node is tested
# against 16129 spherical harmonics, and the symmetry check compares each of the 48 images of
# each node with every node — about 1.5 hours at 100 digits, measured by extrapolation from
# degree 59, where 44% of the time went to the symmetry check alone.
#
# Two facts make most of that unnecessary, without trusting anything the construction did.
#
# Symmetry. Two points of ℝ³ lie in the same O_h orbit exactly when their sorted absolute
# coordinates agree. Sorting the nodes by that key groups them into candidate orbits in
# O(N log N); a group is a complete orbit when its members are distinct and as many as the
# key's orbit size, and the rule is invariant when every group is complete with one weight.
#
# Exactness. For an invariant rule, a polynomial and its group average integrate identically
# both under the rule and over the sphere, so only invariant test functions carry
# information. Group-averaging a monomial xᵅ gives zero when any exponent is odd and depends
# only on the multiset of exponents otherwise, so the test set is the even exponents in
# sorted order. On the sphere, the monomials of one even degree m span all even polynomials
# of degree ≤ m (multiply by (Σxᵢ²)ʲ = 1), and odd polynomials integrate to zero on both
# sides by central symmetry, which O_h contains. So exactness to degree d is a check at the
# single even degree ≤ d, and sharpness one at the next even degree.
#
# At degree 125 that is 352 test functions — the count of the invariants, as it must be —
# plus 363 for sharpness, evaluated once per orbit (132 of them). The test functions are
# monomials evaluated on the delivered nodes, not the p₄, p₆ invariants in orbit parameters
# that the solver used, so the check stays independent of the construction.

"""
    octahedral_orbits(xs, ws, tol) -> Union{Nothing, Tuple}

Group the nodes of a rule on the unit sphere into `O_h` orbits. Returns the orbit
representatives (sorted absolute coordinates, descending) and each orbit's total weight, or
`nothing` if some group is not a complete orbit with a single weight — in which case the
rule is not `O_h`-invariant to tolerance `tol`, or its nodes are too close to tell.
"""
function octahedral_orbits(xs, ws, tol)
    length(first(xs)) == 3 || return nothing
    wscale = maximum(abs, ws)
    keys = [sort!(abs.(collect(x)); rev = true) for x in xs]
    order = sortperm(keys)
    reps = Vector{eltype(keys)}()
    totals = similar(ws, 0)
    i = 1
    while i <= length(order)
        k = keys[order[i]]
        j = i
        while j < length(order) && maximum(abs, keys[order[j + 1]] .- k) <= tol
            j += 1
        end
        grp = order[i:j]
        w = ws[grp[1]]
        all(g -> abs(ws[g] - w) <= tol * wscale, grp) || return nothing
        length(grp) == _octahedral_orbit_size(k, tol) || return nothing
        pts = sort([collect(xs[g]) for g in grp])
        all(m -> maximum(abs, pts[m + 1] .- pts[m]) > tol, 1:(length(pts) - 1)) || return nothing
        push!(reps, k)
        push!(totals, w * length(grp))
        i = j + 1
    end
    return reps, totals
end

"Size of the `O_h` orbit of a point given by its sorted absolute coordinates `k`."
function _octahedral_orbit_size(k, tol)
    zeros_ = count(v -> v <= tol, k)
    same12 = abs(k[1] - k[2]) <= tol
    same23 = abs(k[2] - k[3]) <= tol
    perms = same12 && same23 ? 1 : (same12 || same23) ? 3 : 6
    return perms * 2^(3 - zeros_)
end

"""
    OctahedralTestSet{S}(d, sharp)

The group-averaged even monomials used to verify an `O_h`-symmetric rule of degree `d` on
orbit representatives: all sorted even exponents of the largest even degree `≤ d`, and, for
sharpness, of the next even degree. Evaluated at a representative `k`, a test function is
`(1/6) Σ_σ Πᵢ k_{σ(i)}^{αᵢ}`; with each orbit's total weight this gives exactly the rule's
integral of `xᵅ`.
"""
struct OctahedralTestSet{S} <: VerificationBasis
    exps::Vector{NTuple{3,Int}}
    blockranges::Vector{UnitRange{Int}}
    maxdeg::Int
end
function OctahedralTestSet{S}(d::Integer, sharp::Bool) where {S}
    exps = NTuple{3,Int}[]
    ranges = UnitRange{Int}[]
    me = iseven(d) ? d : d - 1
    for m in (sharp && isodd(d) ? (me, d + 1) : (me,))
        lo = length(exps) + 1
        h = m ÷ 2
        for a in h:-1:0, b in min(a, h - a):-1:0
            c = h - a - b
            0 <= c <= b && push!(exps, (2a, 2b, 2c))
        end
        push!(ranges, lo:length(exps))
    end
    return OctahedralTestSet{S}(exps, ranges, maximum(sum, exps))
end
function evaluate!(b::OctahedralTestSet{S}, k) where {S}
    P = [S(k[j])^e for e in 0:b.maxdeg, j in 1:3]          # P[e+1, j] = k_j^e
    φ = zeros(S, length(b.exps))
    g = zeros(S, length(b.exps))
    for (n, α) in enumerate(b.exps), σ in permutations_of(3)
        t = P[α[1] + 1, σ[1]] * P[α[2] + 1, σ[2]] * P[α[3] + 1, σ[3]]
        φ[n] += t
        # |∇xᵅ| summed componentwise: αᵢ xᵅ / xᵢ, written without dividing by a zero coordinate
        for i in 1:3
            α[i] == 0 && continue
            q = S(α[i])
            for j in 1:3
                q *= P[(j == i ? α[j] - 1 : α[j]) + 1, σ[j]]
            end
            g[n] += q
        end
    end
    return φ ./ 6, g ./ 6
end
blocks(b::OctahedralTestSet) = b.blockranges
exact_integrals(b::OctahedralTestSet{S}) where {S} = S[S(sphere_moment(3, α)) for α in b.exps]

function verify(r::QuadratureRule{D}, c::PolynomialDegree; degree = nothing, bits = nothing,
                integrals = nothing, use_symmetry::Bool = true) where {D}
    d = degree === nothing ? c.d : Int(degree)
    test_sharp = degree === nothing || d == c.d
    exact = _is_exact_type(eltype(r))
    pbits = _rule_precision(r)
    cbits = bits === nothing ? max(2pbits, 128) : Int(bits)
    S = exact ? Rational{BigInt} : BigFloat
    nmax = test_sharp ? d + 1 : d
    return with_bits(cbits) do
        ε = exact ? big(0) : ldexp(big(1.0), -(pbits - 1))
        floor_tol = exact ? big(0.0) : ldexp(big(1.0), -(cbits - 16))
        xs, ws = reference_nodes(r, S)
        refdom = reference(r.domain)
        # an O_h-symmetric rule on S² is checked on its orbit representatives; if its nodes
        # do not group into complete orbits, it falls through to the general check below
        orbits = use_symmetry && !exact && refdom isa Sphere{3} && r.provenance.symmetry === :Oh ?
                 octahedral_orbits(xs, ws, 64ε) : nothing
        if orbits !== nothing
            reps, totals = orbits
            basis = OctahedralTestSet{S}(d, test_sharp)
            bname = "O_h-averaged even monomials on $(length(reps)) orbit representatives (exact sphere moments)"
            res = block_residuals(basis, reps, totals, ε, floor_tol)
            ex = res[1:1]
            sh = test_sharp && isodd(d) ? res[2] : nothing
        elseif refdom isa Simplex && D >= 4
            basis, bname = verification_basis(refdom, test_sharp ? [d, d + 1] : [d], S, exact)
            res = block_residuals(basis, xs, ws, ε, floor_tol)
            ex = res[1:1]
            sh = test_sharp ? res[2] : nothing
        else
            basis, bname = verification_basis(refdom, nmax, S, exact)
            # in one dimension on an orthonormal basis, test sharpness on p_a p_b rather than
            # p_{d+1}; see ProductTop
            if test_sharp && (basis isa RecurrenceBasis || (basis isa JacobiBasis && basis.normalize))
                basis = ProductTop{S}(basis, d + 1)
            end
            res = block_residuals(basis, xs, ws, ε, floor_tol)
            ex = res[1:(d + 1)]
            sh = test_sharp ? res[d + 2] : nothing
        end
        max_res = maximum(t -> t[1], ex)
        max_tol = maximum(t -> t[3], ex)
        is_exact = all(t -> t[2] <= 1, ex)
        sharp = sh === nothing ? nothing : sh[2] > 100 ? true : sh[2] <= 1 ? false : nothing
        # a centrally symmetric rule integrates every odd polynomial exactly, so an even
        # degree claim is understated by construction: exact at d + 1 without testing
        orbits !== nothing && test_sharp && iseven(d) && (sharp = false)
        # structural invariants
        μ = S(measure(refdom))
        wsum = sum(ws)
        wsum_ok = exact ? wsum == μ : abs(wsum - μ) <= 16ε * sum(abs, ws) + floor_tol
        interior = all(x -> isinterior(x, r.domain), _node_iter(r))
        positive = all(>(0), r.weights)
        # the orbit grouping is itself a complete invariance check, done in O(N log N)
        symmetric = orbits !== nothing ? true :
                    check_symmetry(r.provenance.symmetry, xs, ws, exact ? big(0.0) : 64ε)
        Verification(basis = bname, degree = d, max_residual = BigFloat(max_res), tolerance = BigFloat(max_tol),
                     exact = is_exact, sharp = sharp,
                     sharp_residual = sh === nothing ? big(0.0) : BigFloat(sh[1]),
                     weights_sum_ok = wsum_ok, interior = interior, positive = positive, symmetric = symmetric,
                     method = :exact_integration, empirical = false, precision_bits = exact ? 0 : cbits)
    end
end

_node_iter(r::QuadratureRule{1}) = r.nodes
_node_iter(r::QuadratureRule) = r.nodes

"""
    check_simplex_symmetry(xs, ws, tol)

Whether the node/weight set on the reference simplex is invariant under every permutation
of barycentric coordinates, to tolerance `tol` in coordinates and relative weight.
"""
function check_simplex_symmetry(xs, ws, tol)
    N = length(first(xs)) + 1
    perms = permutations_of(N)
    λs = [vcat(1 - sum(x), x) for x in xs]
    wscale = maximum(abs, ws)
    for (λ, w) in zip(λs, ws), σ in perms
        μ = λ[σ]
        found = any(eachindex(λs)) do j
            maximum(abs.(λs[j] .- μ)) <= tol && abs(ws[j] - w) <= tol * wscale
        end
        found || return false
    end
    return true
end

function verify(r::QuadratureRule, c::SpanOf; integrals = nothing, bits = nothing, degree = nothing,
                use_symmetry::Bool = true)
    integrals === nothing && throw(ArgumentError("verifying a SpanOf claim needs the exact `integrals` of the basis"))
    length(integrals) == length(c.basis) || throw(DimensionMismatch("one integral per basis function"))
    exact = _is_exact_type(eltype(r))
    pbits = _rule_precision(r)
    cbits = bits === nothing ? max(2pbits, 128) : Int(bits)
    with_bits(cbits) do
        ε = exact ? big(0) : ldexp(big(1.0), -(pbits - 1))
        mres = big(0.0); mtol = big(0.0); ok = true
        for (b, I) in zip(c.basis, integrals)
            vals = [b(x) for x in r.nodes]
            s = sum(big.(r.weights) .* big.(vals))
            res = abs(s - big(I))
            tol = 16ε * sum(abs.(big.(r.weights)) .* abs.(big.(vals))) + ldexp(big(1.0), -(cbits - 16))
            ok &= res <= tol
            mres = max(mres, res); mtol = max(mtol, tol)
        end
        Verification(basis = "supplied basis ($(length(c.basis)) functions)", degree = -1, max_residual = mres,
                     tolerance = mtol, exact = ok, sharp = nothing, sharp_residual = big(0.0),
                     weights_sum_ok = true, interior = all(x -> isinterior(x, r.domain), r.nodes),
                     positive = all(>(0), r.weights), symmetric = nothing, precision_bits = cbits)
    end
end

"""
    verify_convergence(rules, f, reference; rtol = 0) -> Verification

Empirical check for rules without an exactness claim: the errors `|Q_k f - reference|` over
the sequence `rules` must converge — decreasing until they reach a floor, and never growing
substantially after that — and the last must be within `rtol` of the reference if
`rtol > 0`. The result is flagged `empirical`.

A plateau is accepted, because one is expected: a rule delivered as explicit nodes cannot
resolve an endpoint singularity beyond the point where `1 - x` loses its significant digits,
so tanh-sinh on `1/sqrt(1-x^2)` stalls at about the square root of the working precision.
What the check rules out is divergence.
"""
function verify_convergence(rules::AbstractVector{<:QuadratureRule}, f, reference; rtol = nothing)
    errs = [abs(integrate(f, r) - reference) for r in rules]
    # the default target is the square root of the working precision, the best a rule with
    # explicit nodes can do on an endpoint singularity
    target = rtol === nothing ? 16 * sqrt(eps(float(real(eltype(first(rules)))))) : rtol
    # a step up is allowed only towards the floor: the best error seen, or plain roundoff
    # (a rule can hit an exact zero, which would otherwise make the floor zero)
    ε = eps(float(real(eltype(first(rules)))))
    plateau = max(4 * minimum(errs), 64 * ε * max(abs(reference), one(abs(reference))))
    decreasing = all(k -> errs[k + 1] <= max(errs[k] * 3 // 2, plateau), 1:(length(errs) - 1))
    final_ok = last(errs) <= target * max(abs(reference), one(abs(reference)))
    return Verification(basis = "convergence sweep over $(length(rules)) rules", degree = -1,
                        max_residual = BigFloat(last(errs)),
                        tolerance = BigFloat(target * max(abs(reference), one(abs(reference)))),
                        exact = decreasing && final_ok, sharp = nothing, sharp_residual = big(0.0),
                        weights_sum_ok = true, interior = all(r -> all(x -> isinterior(x, r.domain), r.nodes), rules),
                        positive = all(r -> all(>(0), r.weights), rules), symmetric = nothing,
                        method = :convergence_sweep, empirical = true, precision_bits = 0)
end

"""
    @test_exact rule degree

A `Test.@test` that `rule` integrates every polynomial up to `degree` exactly (to rounding),
for use in users' own test suites.
"""
macro test_exact(r, d)
    ex = :($(check)($(esc(r)); degree = $(esc(d))).exact)
    return Expr(:macrocall, GlobalRef(Test, Symbol("@test")), __source__, ex)
end

function Base.show(io::IO, ::MIME"text/plain", v::Verification)
    status = passed(v) ? "PASSED" : "FAILED"
    println(io, "Verification: ", status, v.empirical ? " (empirical)" : "")
    println(io, "  basis     : ", v.basis, v.precision_bits > 0 ? ", $(v.precision_bits)-bit arithmetic" : ", exact arithmetic")
    v.degree >= 0 && println(io, "  exactness : degree ", v.degree, v.exact ? " ✓" : " ✗",
                              @sprintf("  (max residual %.2e, tolerance %.2e)", v.max_residual, v.tolerance))
    v.sharp === nothing || println(io, "  sharpness : ", v.sharp ? "not exact at degree $(v.degree + 1) ✓" :
                                   "exact at degree $(v.degree + 1) — claim understated ✗",
                                   @sprintf("  (residual %.2e)", v.sharp_residual))
    println(io, "  structure : Σw = measure ", v.weights_sum_ok ? "✓" : "✗",
            ", interior ", v.interior ? "✓" : "✗", ", positive ", v.positive ? "✓" : "✗",
            v.symmetric === nothing ? "" : ", symmetric " * (v.symmetric ? "✓" : "✗"))
end

"""
    check_symmetry(group, xs, ws, tol)

Whether the rule is invariant under the symmetry group it claims: `:S_N` permutes the
barycentric coordinates of a simplex, `:reflection` is `x ↦ -x` on the reference interval,
`:Oh` is the 48 signed permutations of Cartesian coordinates on a sphere. `nothing` when no
symmetry is claimed.
"""
check_symmetry(group::Symbol, xs, ws, tol) =
    group === :none ? nothing :
    group === :reflection ? check_reflection_symmetry(xs, ws, tol) :
    group === :Oh ? check_octahedral_symmetry(xs, ws, tol) :
    check_simplex_symmetry(xs, ws, tol)

"""
    check_octahedral_symmetry(xs, ws, tol)

Whether the node/weight set is invariant under the full octahedral group: every image of
every node under a signed permutation of its coordinates is again a node, carrying the same
weight. Checking the 48 images of each node is what makes this independent of the orbit
machinery that built the rule.
"""
function check_octahedral_symmetry(xs, ws, tol)
    length(first(xs)) == 3 || return false
    wscale = maximum(abs, ws)
    for (x, w) in zip(xs, ws), p in permutations_of(3), s1 in (1, -1), s2 in (1, -1), s3 in (1, -1)
        y = (s1 * x[p[1]], s2 * x[p[2]], s3 * x[p[3]])
        any(eachindex(xs)) do j
            maximum(abs, (xs[j][1] - y[1], xs[j][2] - y[2], xs[j][3] - y[3])) <= tol &&
                abs(ws[j] - w) <= tol * wscale
        end || return false
    end
    return true
end

"Whether the node/weight set on the reference interval is invariant under `x ↦ -x`."
function check_reflection_symmetry(xs, ws, tol)
    wscale = maximum(abs, ws)
    for (x, w) in zip(xs, ws)
        any(j -> abs(xs[j][1] + x[1]) <= tol && abs(ws[j] - w) <= tol * wscale, eachindex(xs)) || return false
    end
    return true
end
