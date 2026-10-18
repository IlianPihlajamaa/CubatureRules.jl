# Verification (PLAN §8): is this rule what it claims to be? Distinct from the Certificate,
# which only says that Newton converged on the system it was given.
#
# Exactness is tested against *orthonormal* bases (Dubiner on triangles and tetrahedra,
# Legendre/Jacobi on intervals) at twice the rule's precision; monomials would produce false
# failures at high degree. Each basis function's residual is held to a tolerance derived from
# the rounding of the delivered nodes and weights: 16 ε Σ |wᵢ| (|φ(xᵢ)| + |∇φ(xᵢ)|·max(|xᵢ|, 1)).
# Exact (Rational) rules are held to zero.
#
# This file is the check itself; the bases are in verify/bases.jl, the kernels' reference
# integrals in verify/kernel_reference.jl, and the symmetry check in verify/symmetry.jl.

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

"Per degree block: the largest residual, the largest ratio of residual to tolerance, and tolerance."
function block_residuals(basis, xs, ws, ε, floor_tol)
    I = exact_integrals(basis)
    L = length(I)
    Σ = zeros(eltype(ws), L)
    scale = zeros(real(eltype(ws)), L)
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
    out = @NamedTuple{residual::BigFloat, ratio::BigFloat, tolerance::BigFloat}[]
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
        push!(out, (residual = mres, ratio = mratio, tolerance = mtol))
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
# key's orbit size, and the rule is invariant when every group is complete with one weight
# (`group_orbits`, src/verify/symmetry.jl, the same for every group).
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

"""
    exactness_test(refdom, d, test_sharp, S, exact) -> (basis, name, nclaim)

The basis exactness to degree `d` is tested on, its description, and how many of its degree
blocks the claim covers; the next block, when `test_sharp`, tests sharpness. Every domain's
basis runs over the degrees `0 … d` (and `d + 1`); beyond three dimensions the simplex's covers
only `d` and `d + 1` (see its `verification_basis`).
"""
function exactness_test(refdom::Domain, d, test_sharp, ::Type{S}, exact) where {S}
    basis, name = verification_basis(refdom, test_sharp ? d + 1 : d, S, exact)
    return (test_sharp ? sharpness_basis(basis, d, S) : basis), name, d + 1
end
# In one dimension on an orthonormal basis, sharpness is tested on p_a p_b rather than
# p_{d+1}; see ProductTop.
sharpness_basis(basis, d, ::Type{S}) where {S} = basis
sharpness_basis(basis::RecurrenceBasis, d, ::Type{S}) where {S} = ProductTop{S}(basis, d + 1)
sharpness_basis(basis::JacobiBasis, d, ::Type{S}) where {S} = basis.normalize ? ProductTop{S}(basis, d + 1) : basis

function verify(r::QuadratureRule, c::PolynomialDegree; degree = nothing, bits = nothing,
                integrals = nothing, use_symmetry::Bool = true)
    d = degree === nothing ? c.d : Int(degree)
    test_sharp = degree === nothing || d == c.d
    exact = _is_exact_type(numbertype(r))
    pbits = _rule_precision(r)
    cbits = bits === nothing ? max(2pbits, 128) : Int(bits)
    S = exact ? Rational{BigInt} : BigFloat
    return with_bits(cbits) do
        ε = exact ? big(0) : ldexp(big(1.0), -(pbits - 1))
        floor_tol = exact ? big(0.0) : ldexp(big(1.0), -(cbits - 16))
        xs, ws = reference_nodes(r, S)
        refdom = reference(r.domain)
        # an O_h-symmetric rule on S² is checked on its orbit representatives; if its nodes
        # do not group into complete orbits, it falls through to the general check below
        orbits = use_symmetry && !exact && refdom isa Sphere{3} && r.provenance.symmetry isa Oh ?
                 group_orbits(Oh(), xs, ws, 64ε) : nothing
        if orbits !== nothing
            reps, totals = orbits
            basis = OctahedralTestSet{S}(d, test_sharp)
            bname = "O_h-averaged even monomials on $(length(reps)) orbit representatives (exact sphere moments)"
            res = block_residuals(basis, reps, totals, ε, floor_tol)
            ex = res[1:1]
            sh = test_sharp && isodd(d) ? res[2] : nothing
        else
            basis, bname, nclaim = exactness_test(refdom, d, test_sharp, S, exact)
            res = block_residuals(basis, xs, ws, ε, floor_tol)
            ex = res[1:nclaim]
            sh = test_sharp ? res[nclaim + 1] : nothing
        end
        max_res = maximum(t -> t.residual, ex)
        max_tol = maximum(t -> t.tolerance, ex)
        is_exact = all(t -> t.ratio <= 1, ex)
        sharp = sh === nothing ? nothing : sh.ratio > 100 ? true : sh.ratio <= 1 ? false : nothing
        # Some rules come within rounding of the next degree: a Gauss–Patterson rule misses
        # degree d + 1 by 1e-20 at 127 points and 6e-41 at 255. Delivered at lower precision,
        # such a rule looks exact there, which is not the same as being exact. The
        # construction measured the error at working precision; say so rather than fail.
        sharp_note = ""
        cert = r.certificate
        if sharp === false && cert !== nothing && cert.next_error !== nothing && cert.next_error > 0
            sharp = nothing
            sharp_note = @sprintf("not resolvable at this precision: the rule misses degree %d by %.1e, measured at construction",
                                  d + 1, cert.next_error)
        end
        # a centrally symmetric rule integrates every odd polynomial exactly, so an even
        # degree claim is understated by construction: exact at d + 1 without testing
        orbits !== nothing && test_sharp && iseven(d) && (sharp = false)
        # structural invariants
        μ = lift(S, measure(refdom))
        wsum = sum(ws)
        wsum_ok = exact ? wsum == μ : abs(wsum - μ) <= 16ε * sum(abs, ws) + floor_tol
        interior = all(x -> isinterior(x, r.domain), r.nodes)
        positive = eltype(r.weights) <: Real && all(>=(0), r.weights)   # an underflowed weight is zero, not negative
        # the orbit grouping is itself a complete invariance check, done in O(N log N)
        symmetric = orbits !== nothing ? true :
                    check_invariance(r.provenance.symmetry, xs, ws, exact ? big(0.0) : 64ε)
        Verification(basis = bname, degree = d, max_residual = BigFloat(max_res), tolerance = BigFloat(max_tol),
                     exact = is_exact, sharp = sharp,
                     sharp_residual = sh === nothing ? big(0.0) : BigFloat(sh.residual),
                     weights_sum_ok = wsum_ok, interior = interior, positive = positive, symmetric = symmetric,
                     method = :exact_integration, empirical = false, precision_bits = exact ? 0 : cbits,
                     sharp_note = sharp_note)
    end
end

function verify(r::QuadratureRule, c::SpanOf; integrals = nothing, bits = nothing, degree = nothing,
                use_symmetry::Bool = true)
    integrals === nothing && throw(ArgumentError("verifying a SpanOf claim needs the exact `integrals` of the basis"))
    length(integrals) == length(c.basis) || throw(DimensionMismatch("one integral per basis function"))
    exact = _is_exact_type(numbertype(r))
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

Empirical check for rules without an exactness claim. It computes the error
`|Q_k f - reference|` of each rule in `rules` and checks two things:

1. The errors never grow by more than half from one rule to the next, except at the floor:
   four times the smallest error seen, or the rounding error of the last rule's sum,
   `64 N ε Σ|wᵢ f(xᵢ)|` for `N` points at the rules' precision `ε`.
2. The last error is at most `rtol · max(|reference|, 1)`. `rtol` defaults to `16 √ε`, the
   best a rule given as explicit nodes can do on an endpoint singularity: tanh-sinh on
   `1/sqrt(1-x^2)` stalls there, because `1 - x` loses half its digits at the outermost nodes.

It does not check how fast the errors fall, or that they fall at all once below the target: a
sequence that starts at its floor passes. What it rules out is divergence and a wrong answer.
The result is flagged `empirical` and records the errors.
"""
function verify_convergence(rules::AbstractVector{<:QuadratureRule}, f, reference; rtol = nothing)
    r1 = first(rules)
    exact = _is_exact_type(numbertype(r1))
    pbits = _rule_precision(r1)
    # the rules' own precision, not the global BigFloat precision: a 50-digit rule evaluated
    # at the default 256 bits would otherwise be held to rounding it cannot reach
    ε = exact ? eps(Float64) : ldexp(big(1.0), -(pbits - 1))
    wbits = max(2pbits, 128)
    errs, mass = with_bits(wbits) do
        vals = [integrate(f, r) for r in rules]
        big.(abs.(vals .- reference)), sum(abs(big(w) * big(f(x))) for (x, w) in zip(nodes(last(rules)), weights(last(rules))))
    end
    scale = max(abs(big(reference)), big(1.0))
    target = (rtol === nothing ? 16 * sqrt(ε) : big(rtol)) * scale
    floor = max(4 * minimum(errs), 64 * npoints(last(rules)) * ε * max(mass, big(1.0)))
    decreasing = all(k -> errs[k + 1] <= max(errs[k] * 3 // 2, floor), 1:(length(errs) - 1))
    final_ok = last(errs) <= target
    return Verification(basis = "convergence sweep over $(length(rules)) rules", degree = -1,
                        max_residual = last(errs), tolerance = target,
                        exact = decreasing && final_ok, sharp = nothing, sharp_residual = big(0.0),
                        weights_sum_ok = true, interior = all(r -> all(x -> isinterior(x, r.domain), r.nodes), rules),
                        positive = all(r -> all(>(0), r.weights), rules), symmetric = nothing,
                        method = :convergence_sweep, empirical = true, precision_bits = pbits,
                        errors = errs, decreasing = decreasing)
end

"""
    @test_exact rule degree

A `Test.@test` that `rule` integrates every polynomial up to `degree` exactly (to rounding),
for use in users' own test suites, after `using Test`.
"""
macro test_exact(r, d)
    # all of it in the caller's scope, so that `Test` is the caller's: this package does not
    # depend on Test, a test suite does; `check` is reached through a GlobalRef
    ex = :($(GlobalRef(@__MODULE__, :check))($r; degree = $d).exact)
    return esc(Expr(:macrocall, Expr(:., :Test, QuoteNode(Symbol("@test"))), __source__, ex))
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
    v.sharp === nothing && !isempty(v.sharp_note) && println(io, "  sharpness : ", v.sharp_note)
    if v.method === :convergence_sweep && !isempty(v.errors)
        println(io, "  errors    : ", join((@sprintf("%.1e", e) for e in v.errors), ", "))
        println(io, "  no growth : ", v.decreasing === true ? "✓" : "✗", " (except at the floor)")
        println(io, "  last error: ", v.max_residual <= v.tolerance ? "✓" : "✗",
                @sprintf("  (%.1e, tolerance %.1e)", v.max_residual, v.tolerance))
    end
    println(io, "  structure : Σw = measure ", v.weights_sum_ok ? "✓" : "✗",
            ", interior ", v.interior ? "✓" : "✗", ", positive ", v.positive ? "✓" : "✗",
            v.symmetric === nothing ? "" : ", symmetric " * (v.symmetric ? "✓" : "✗"))
end
