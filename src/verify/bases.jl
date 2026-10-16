# The verification bases (PLAN §8): orthonormal polynomials on each reference domain, the
# exact integrals they are checked against, and the nodes of a rule mapped to the reference
# domain they are evaluated at. The kernels singular at a point of a triangle or a
# tetrahedron have their own, with integrals computed in polar coordinates
# (verify/kernel_reference.jl); the checks themselves are in verify/verify.jl.

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
    x = S(xv[1])
    p, dp = zeros(S, b.n + 1), zeros(S, b.n + 1)
    if b.normalize
        orthonormal!(p, dp, jacobi_recurrence(b.α, b.β), x)
    else   # exact arithmetic: unnormalised Legendre, α = β = 0
        three_term!(p, dp, x, one(S), k -> (2k + 1, 0, k, k + 1))
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

# The reference wedge: the orthonormal Dubiner polynomials of the triangle times orthonormal
# Legendre polynomials in z, graded by total degree. Orthonormal on the wedge, so only the
# constant has a nonzero integral.
struct WedgeBasis{S} <: VerificationBasis
    tri::DubinerBasis{S}
    line::JacobiBasis{S,Int}
    idx::Vector{Tuple{Int,Int}}                 # (Dubiner index, Legendre degree)
    blockranges::Vector{UnitRange{Int}}
end
function WedgeBasis{S}(n::Int) where {S}
    idx, ranges = Tuple{Int,Int}[], UnitRange{Int}[]
    for k in 0:n
        lo = length(idx) + 1
        for r in 0:k, p in 0:(k - r)
            push!(idx, (dubiner_index(p, k - r - p), r))
        end
        push!(ranges, lo:length(idx))
    end
    return WedgeBasis{S}(DubinerBasis{S}(n; normalize = true), JacobiBasis{S,Int}(n, 0, 0, true), idx, ranges)
end
function evaluate!(b::WedgeBasis{S}, x) where {S}
    φt, gt = evaluate!(b.tri, (x[1], x[2]))
    φl, gl = evaluate!(b.line, (x[3],))
    return [φt[i] * φl[r + 1] for (i, r) in b.idx], [gt[i] * abs(φl[r + 1]) + abs(φt[i]) * gl[r + 1] for (i, r) in b.idx]
end
blocks(b::WedgeBasis) = b.blockranges
function exact_integrals(b::WedgeBasis{S}) where {S}
    v = zeros(S, length(b.idx))
    v[1] = exact_integrals(b.tri)[1] * exact_integrals(b.line)[1]
    return v
end

# The reference pyramid: orthonormal Legendre polynomials on its bounding box
# [-1, 1]² × [0, 1], graded by total degree. They are not orthogonal on the pyramid but are
# of order one on it, which monomials are not. Their integrals over the pyramid use
#
#     ∫_{-s}^{s} P_a = 2 (P_{a+1}(s) − P_{a−1}(s)) / (2a + 1) for even a ≥ 2, 2s for a = 0,
#
# zero for odd a, on each slice z = 1 − s, and a Gauss–Legendre rule in z exact for the
# degree. The conical rules collapse the pyramid onto a cube instead, so the two share
# nothing but one-dimensional Gauss nodes.
struct PyramidBasis{S} <: VerificationBasis
    n::Int
    axis::JacobiBasis{S,Int}
    idx::Vector{NTuple{3,Int}}
    blockranges::Vector{UnitRange{Int}}
end
function PyramidBasis{S}(n::Int) where {S}
    idx, ranges = NTuple{3,Int}[], UnitRange{Int}[]
    for k in 0:n
        lo = length(idx) + 1
        for c in compositions(k, 3)
            push!(idx, NTuple{3,Int}(c))
        end
        push!(ranges, lo:length(idx))
    end
    return PyramidBasis{S}(n, JacobiBasis{S,Int}(n, 0, 0, true), idx, ranges)
end
function evaluate!(b::PyramidBasis{S}, x) where {S}
    px, gx = evaluate!(b.axis, (x[1],))
    py, gy = evaluate!(b.axis, (x[2],))
    pz, gz = evaluate!(b.axis, (2x[3] - 1,))
    φ = [px[i + 1] * py[j + 1] * pz[k + 1] for (i, j, k) in b.idx]
    g = [gx[i + 1] * abs(py[j + 1] * pz[k + 1]) + gy[j + 1] * abs(px[i + 1] * pz[k + 1]) +
         2gz[k + 1] * abs(px[i + 1] * py[j + 1]) for (i, j, k) in b.idx]
    return φ, g
end
blocks(b::PyramidBasis) = b.blockranges
function exact_integrals(b::PyramidBasis{S}) where {S}
    n = b.n
    M = cld(n + 3, 2)                            # the integrand in z has degree ≤ n + 2
    t, ω, _ = gauss_jacobi_work(M, 0, 0, precision(S))
    nrm(a) = sqrt(S(2a + 1) / 2)                 # orthonormal Legendre is nrm(a) P_a
    F = zeros(S, n + 1, M)                       # F[a + 1, m] = ∫_{-s}^{s} p_a, s = 1 − z_m
    pz = zeros(S, n + 1, M)
    for m in 1:M
        z = (1 + S(t[m])) / 2
        s = 1 - z
        P = zeros(S, n + 2)                      # P[a + 1] = P_a(s)
        P[1] = one(S)
        n >= 0 && (P[2] = s)
        for a in 1:n
            P[a + 2] = ((2a + 1) * s * P[a + 1] - a * P[a]) / (a + 1)
        end
        for a in 0:n
            F[a + 1, m] = a == 0 ? 2s * nrm(0) : isodd(a) ? zero(S) : 2 * (P[a + 2] - P[a]) / (2a + 1) * nrm(a)
        end
        pz[:, m] = evaluate!(b.axis, (S(t[m]),))[1]
    end
    return [sum(S(ω[m]) / 2 * F[i + 1, m] * F[j + 1, m] * pz[k + 1, m] for m in 1:M) for (i, j, k) in b.idx]
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
    p, dp = orthonormal!(zeros(S, b.n + 1), zeros(S, b.n + 1), b.rec, S(xv[1]))
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

# A wedge or a pyramid through its affine frame
function reference_nodes(r::QuadratureRule{3,T,<:Union{Wedge,Pyramid}}, ::Type{S}) where {T,S}
    A, b = affine_frame(r.domain, S)
    J = abs(det(A))
    xs = [Vector{S}(A \ (SVector{3,S}(map(S, x)) - b)) for x in r.nodes]
    return xs, [S(w) / J for w in r.weights]
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

# A triangle in space maps back to the reference triangle through the left inverse of its
# frame, and its weights carry the area element |e₁ × e₂|.
function reference_nodes(r::QuadratureRule{3,T,<:SurfaceTriangle}, ::Type{S}) where {T,S}
    E, v0, J = surface_frame(r.domain, S)
    P = inv(E' * E) * E'
    xs = [Vector{S}(P * (SVector{3,S}(map(S, x)) - v0)) for x in r.nodes]
    return xs, [S(w) / J for w in r.weights]
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
# and so is a weight given as a function, and a singular kernel
reference_nodes(r::QuadratureRule{1,T,<:FunctionDomain}, ::Type{S}) where {T,S} =
    ([[S(x)] for x in r.nodes], [S(w) for w in r.weights])
reference_nodes(r::QuadratureRule{1,T,<:KernelDomain}, ::Type{S}) where {T,S} =
    ([[S(x)] for x in r.nodes], [S(w) for w in r.weights])
# and an oscillatory weight, whose weights may be complex
reference_nodes(r::QuadratureRule{1,T,<:OscillatoryDomain}, ::Type{S}) where {T,S} =
    ([[S(x)] for x in r.nodes], [lift(S, w) for w in r.weights])

"A real or complex number in the verification type `S`, or `Complex{S}`."
lift(::Type{S}, z::Complex) where {S} = Complex{S}(z)
lift(::Type{S}, x) where {S} = S(x)

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
verification_basis(dom::Wedge, n, ::Type{S}, exact) where {S} = WedgeBasis{S}(n),
    "orthonormal Dubiner on the triangle × orthonormal Legendre in z, graded by total degree"
verification_basis(dom::Pyramid, n, ::Type{S}, exact) where {S} = PyramidBasis{S}(n),
    "orthonormal Legendre on the bounding box, graded by total degree (integrals over the pyramid by slices)"
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

verification_basis(dom::ExponentialDomain{D}, n, ::Type{S}, exact) where {D,S} =
    MonomialTestSet{D,S}(_degrees(n), α -> exponential_moment(D, α)),
    "monomials of degree $(join(n, ", ")) (exact e^-|x| moments)"

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
    N = maximum(b.degs) + 1
    aux = b.w.aux                                    # monic: α = d = 1, c_k = β_k
    p, dp = three_term!(zeros(S, N), zeros(S, N), S(xv[1]), one(S),
                        k -> (one(S), S(aux.a(k, S)), S(aux.b(k, S)), one(S)))
    return S[p[k + 1] for k in b.degs], S[abs(dp[k + 1]) for k in b.degs]
end
blocks(b::MomentBasis) = [k:k for k in eachindex(b.degs)]
exact_integrals(b::MomentBasis{S}) where {S} = S[S(b.w.moments(k, S)) for k in b.degs]

verification_basis(dom::MomentDomain, n, ::Type{S}, exact) where {S} =
    MomentBasis{S}(_degrees(n), dom.weight),
    "monic auxiliary polynomials of $(dom.weight.label), against the given moments"

# A weight given as a function is verified against its modified moments, each computed by
# discretising the weight until two resolutions agree at the verification precision. That
# shares the discretisation with the construction but not the Stieltjes procedure, the
# recurrence or the Newton iteration.
verification_basis(dom::FunctionDomain, n, ::Type{S}, exact) where {S} =
    MomentBasis{S}(_degrees(n), verification_weight(dom.weight)),
    "monic polynomials on the support, against moments of the weight computed by discretisation"

# The orthogonal polynomials of a Jacobi weight on an interval c ± h, scaled to q₀ = 1 so
# that they stay of order one at every degree (monic ones on [0, 1] fall below any absolute
# tolerance by degree 100), with their integrals against whatever the domain carries.
struct IntervalMomentBasis{S,M} <: VerificationBasis
    degs::Vector{Int}
    a::Vector{S}                  # orthonormal recurrence on [-1, 1]; b[j + 1] = b_j
    b::Vector{S}
    c::S                          # the interval is c ± h
    h::S
    moments::Vector{M}            # S, or Complex{S} for a complex weight
end
function evaluate!(B::IntervalMomentBasis{S}, xv) where {S}
    x = (S(xv[1]) - B.c) / B.h
    N = maximum(B.degs) + 1
    p, dp = three_term!(zeros(S, N), zeros(S, N), x, one(S), k -> (one(S), B.a[k + 1], B.b[k + 1], B.b[k + 2]))
    return S[p[j + 1] for j in B.degs], S[abs(dp[j + 1]) / B.h for j in B.degs]
end
blocks(B::IntervalMomentBasis) = [j:j for j in eachindex(B.degs)]
exact_integrals(B::IntervalMomentBasis{S,M}) where {S,M} = M[B.moments[j + 1] for j in B.degs]

# A principal value or finite part is verified against the kernel applied to these
# polynomials, each computed by Gauss–Jacobi quadrature of a divided difference (see
# `kernel_moments`). Only the Hilbert transform of the weight is shared with the
# construction, which uses a recurrence instead.
function verification_basis(dom::KernelDomain, n, ::Type{S}, exact) where {S}
    degs = collect(_degrees(n))
    K = maximum(degs)
    k = dom.weight
    rec = jacobi_recurrence(k.α, k.β)
    lo, hi = S(dom.base.a), S(dom.base.b)
    basis = IntervalMomentBasis{S,S}(degs, S[rec.a(j, S) for j in 0:K], S[j == 0 ? zero(S) : rec.b(j, S) for j in 0:(K + 1)],
                           (hi + lo) / 2, (hi - lo) / 2, S.(kernel_moments(dom, K)))
    return basis, "orthogonal Jacobi polynomials on the interval, against the kernel's integrals of them " *
                  "(Gauss–Jacobi quadrature of divided differences)"
end

# An oscillatory weight is verified against its Legendre moments, with the spherical Bessel
# functions computed downwards by Miller's algorithm (`oscillatory_moments`), where the
# construction computes them upwards. q_k = √(2k + 1) P_k in the orthonormal scaling.
function verification_basis(dom::OscillatoryDomain, n, ::Type{S}, exact) where {S}
    degs = collect(_degrees(n))
    K = maximum(degs)
    rec = jacobi_recurrence(0, 0)
    lo, hi = S(dom.base.a), S(dom.base.b)
    m = [sqrt(S(2k + 1)) * lift(S, mk) for (k, mk) in zip(0:K, oscillatory_moments(dom, K))]
    basis = IntervalMomentBasis{S,eltype(m)}(degs, S[rec.a(j, S) for j in 0:K],
                                             S[j == 0 ? zero(S) : rec.b(j, S) for j in 0:(K + 1)],
                                             (hi + lo) / 2, (hi - lo) / 2, m)
    return basis, "orthogonal Legendre polynomials on the interval, against the weight's moments " *
                  "(spherical Bessel functions by Miller's algorithm)"
end
