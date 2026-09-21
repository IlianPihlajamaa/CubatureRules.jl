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
"""
verify(r::QuadratureRule; kw...) = verify(r, r.exactness; kw...)

"""
    check(rule; degree, bits) -> Verification

Re-run verification on demand, optionally at a user-specified degree and precision. The
public entry point for rules constructed by hand or loaded from elsewhere.
"""
check(r::QuadratureRule; degree = nothing, bits = nothing) = verify(r; degree, bits)

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

function reference_nodes(r::QuadratureRule{1,T,<:WeightedDomain}, ::Type{S}) where {T,S}
    isreference(r.domain.base) || throw(ArgumentError("verification of weighted rules needs the reference interval"))
    return [[S(x)] for x in r.nodes], [S(w) for w in r.weights]
end

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
verification_basis(dom::Interval, n, ::Type{S}, exact) where {S} =
    JacobiBasis{S,Int}(n, 0, 0, !exact), exact ? "Legendre (unnormalised, exact arithmetic)" : "orthonormal Legendre"
verification_basis(dom::WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight}, n, ::Type{S}, exact) where {S} =
    JacobiBasis{S,typeof(dom.weight.α)}(n, dom.weight.α, dom.weight.β, true),
    "orthonormal Jacobi($(dom.weight.α), $(dom.weight.β))"

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

function verify(r::QuadratureRule{D}, c::PolynomialDegree; degree = nothing, bits = nothing,
                integrals = nothing) where {D}
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
        if refdom isa Simplex && D >= 4
            basis, bname = verification_basis(refdom, test_sharp ? [d, d + 1] : [d], S, exact)
            res = block_residuals(basis, xs, ws, ε, floor_tol)
            ex = res[1:1]
            sh = test_sharp ? res[2] : nothing
        else
            basis, bname = verification_basis(refdom, nmax, S, exact)
            res = block_residuals(basis, xs, ws, ε, floor_tol)
            ex = res[1:(d + 1)]
            sh = test_sharp ? res[d + 2] : nothing
        end
        max_res = maximum(t -> t[1], ex)
        max_tol = maximum(t -> t[3], ex)
        is_exact = all(t -> t[2] <= 1, ex)
        sharp = sh === nothing ? nothing : sh[2] > 100 ? true : sh[2] <= 1 ? false : nothing
        # structural invariants
        μ = S(measure(refdom))
        wsum = sum(ws)
        wsum_ok = exact ? wsum == μ : abs(wsum - μ) <= 16ε * sum(abs, ws) + floor_tol
        interior = all(x -> isinterior(x, r.domain), _node_iter(r))
        positive = all(>(0), r.weights)
        symmetric = check_symmetry(r.provenance.symmetry, xs, ws, exact ? big(0.0) : 64ε)
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

function verify(r::QuadratureRule, c::SpanOf; integrals = nothing, bits = nothing, degree = nothing)
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

Empirical check for rules without an exactness claim: the errors `|Q_k f - reference|`
over the sequence `rules` must decrease (non-strictly, allowing for rounding at the end)
and the last must be within `rtol` of the reference if `rtol > 0`. The result is flagged
`empirical`.
"""
function verify_convergence(rules::AbstractVector{<:QuadratureRule}, f, reference; rtol = 0)
    errs = [abs(integrate(f, r) - reference) for r in rules]
    floor = 64 * eps(float(real(eltype(first(rules))))) * abs(reference)
    decreasing = all(k -> errs[k + 1] <= errs[k] || errs[k + 1] <= floor, 1:(length(errs) - 1))
    final_ok = rtol == 0 || last(errs) <= rtol * abs(reference)
    return Verification(basis = "convergence sweep over $(length(rules)) rules", degree = -1,
                        max_residual = BigFloat(last(errs)), tolerance = BigFloat(rtol * abs(reference)),
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
barycentric coordinates of a simplex, `:reflection` is `x ↦ -x` on the reference interval.
`nothing` when no symmetry is claimed.
"""
check_symmetry(group::Symbol, xs, ws, tol) =
    group === :none ? nothing :
    group === :reflection ? check_reflection_symmetry(xs, ws, tol) :
    check_simplex_symmetry(xs, ws, tol)

"Whether the node/weight set on the reference interval is invariant under `x ↦ -x`."
function check_reflection_symmetry(xs, ws, tol)
    wscale = maximum(abs, ws)
    for (x, w) in zip(xs, ws)
        any(j -> abs(xs[j][1] + x[1]) <= tol && abs(ws[j] - w) <= tol * wscale, eachindex(xs)) || return false
    end
    return true
end
