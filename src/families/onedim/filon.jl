# Filon-type rules for oscillatory weights (PLAN §6 Tier 2, v0.5). The weights and their
# moments are in domains/oscillatory.jl.
#
# The rule is interpolatory on the N + 1 Gauss–Lobatto points: the end points and the zeros
# of P_N'. With the end points among the nodes, the error of an interpolatory rule for
# e^{iωx} falls like ω⁻² as ω grows at fixed N (Iserles & Nørsett 2005), and as N grows it
# converges like the interpolant, so the same rule serves small and large ω. The Lobatto
# points also make the weights cheap: the Lobatto rule integrates P_j P_k exactly except
# for j = k = N, where it gives 2/N instead of 2/(2N + 1), so the interpolant's Legendre
# coefficients are discrete inner products and
#
#     w_j = λ_j Σ_k P_k(x_j) M_k / γ_k,    γ_k = 2/(2k + 1) for k < N, γ_N = 2/N,
#
# with λ_j the Lobatto weights and M_k = 2 iᵏ j_k(κ) the moments. No linear system is
# solved, and nothing cancels as ω → 0: the moments of high degree are small there, not the
# differences of large numbers that the classical closed-form Filon weights are.

"""
    Filon()

Filon-type rules for an [`Oscillatory`](@ref) weight `e^{iωx}`, `cos(ωx)` or `sin(ωx)` on an
interval: interpolatory on the `N + 1` Gauss–Lobatto points, exact to degree `N` for every
`ω`. Including the end points makes the error fall like `ω⁻²` as `ω` grows at fixed `N`.
The weights come from the Legendre moments of the weight, `2 iᵏ j_k(ωh)` with `j_k` the
spherical Bessel functions, so small `ω` loses no digits to cancellation.

The node set includes both end points, so these rules are not interior.
"""
struct Filon <: RuleFamily end

derivation(::Type{Filon}) = Derived()
family_name(::Filon) = "Filon"

const FILON_1928 = Citation(key = "Filon1928", authors = ["L. N. G. Filon"],
                            title = "On a quadrature formula for trigonometric integrals",
                            journal = "Proceedings of the Royal Society of Edinburgh", year = 1928, volume = "49",
                            pages = "38--47", doi = "10.1017/S0370164600026262")
const ISERLES_NORSETT_2005 = Citation(key = "IserlesNorsett2005", authors = ["Arieh Iserles", "Syvert P. Nørsett"],
                                      title = "Efficient quadrature of highly oscillatory integrals using derivatives",
                                      journal = "Proceedings of the Royal Society A", year = 2005, volume = "461",
                                      pages = "1383--1399", doi = "10.1098/rspa.2004.1401")

candidates(::Type{Filon}, dom::OscillatoryDomain, ::PolynomialDegree) = [Filon()]
filon_order(degree) = max(degree, 1)
npoints(::Filon, dom, degree::Integer) = filon_order(degree) + 1
claimed_degree(::Filon, dom, degree) = filon_order(degree)
degree_range(::Filon, dom) = 0:typemax(Int)
degree_for_npoints(::Filon, dom, n::Integer) = n - 1
properties(::Filon, dom, degree) = (positive = false, interior = false, symmetry = :none, nested = false)

# Legendre values P[k + 1, j] = P_k(x_j), k = 0 … K
function legendre_table(K, x)
    P = Matrix{BigFloat}(undef, K + 1, length(x))
    for (j, xj) in enumerate(x)
        P[1, j] = one(xj)
        K >= 1 && (P[2, j] = xj)
        for k in 1:(K - 1)
            P[k + 2, j] = ((2k + 1) * xj * P[k + 1, j] - k * P[k, j]) / (k + 1)
        end
    end
    return P
end

function build(f::Filon, dom::OscillatoryDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("Filon weights involve sin and cos; $(T) is not supported"))
    check_oscillatory(dom)
    w = dom.weight
    N = filon_order(degree)
    wbits = ctx.bits + gj_guard_bits(N) + 16
    with_bits(wbits) do
        h, c, κ = oscillatory_frame(dom)
        xi = N >= 2 ? first(gauss_jacobi_work(N - 1, 1, 1, wbits)) : BigFloat[]
        x = vcat(-one(BigFloat), BigFloat.(xi), one(BigFloat))
        K = N + 2                                        # one past any claim, to measure the miss
        P = legendre_table(K, x)
        λ = [2 / (N * (N + 1) * P[N + 1, j]^2) for j in eachindex(x)]
        jb = spherical_bessel_up(K, κ)
        M = [2 * im^k * jb[k + 1] for k in 0:K]          # ∫ P_k e^{iκx} on [−1, 1]
        γ = [k < N ? BigFloat(2) / (2k + 1) : BigFloat(2) / N for k in 0:N]
        phase = h * cis(BigFloat(w.ω) * c)
        W = [oscillatory_part(w.kind, phase * λ[j] * sum(P[k + 1, j] * M[k + 1] / γ[k + 1] for k in 0:N))
             for j in eachindex(x)]
        m = [oscillatory_part(w.kind, phase * Mk) for Mk in M]
        # Exact to degree N. A real kernel with a parity on an interval centred at 0 gives one
        # more: cos is even, so every odd polynomial integrates to zero, and the symmetric
        # nodes with symmetric weights give zero too; sin likewise for even polynomials.
        d = N
        iszero(c) && ((w.kind === :cos && iseven(N)) || (w.kind === :sin && isodd(N))) && (d += 1)
        # check the claim at working precision, and measure the miss one degree above it
        tol = ldexp(BigFloat(1), -(wbits - 24))
        miss = zeros(BigFloat, K + 1)
        for k in 0:K
            s = sum(W[j] * P[k + 1, j] for j in eachindex(x))
            scale = sum(abs(W[j] * P[k + 1, j]) for j in eachindex(x)) + abs(m[k + 1])
            miss[k + 1] = abs(s - m[k + 1])
            k <= d && miss[k + 1] > tol * scale &&
                throw(RefinementError("Filon", "the rule for $(w) misses degree $k by $(Float64(miss[k + 1]))"))
        end
        fin(z) = z isa Complex ? complex(finalize_number(ctx, real(z)), finalize_number(ctx, imag(z))) :
                 finalize_number(ctx, z)
        xs = [finalize_number(ctx, c + h * xj) for xj in x]
        ws = [fin(Wj) for Wj in W]
        # the defining equations at the delivered interior nodes, back on [−1, 1]
        res = N >= 2 ? gj_residual([(BigFloat(xk) - c) / h for xk in xs[2:(end - 1)]], 1, 1, wbits) : zero(h)
        # the miss at d + 1 in the basis verification uses, √(2k + 1) P_k
        next = sqrt(BigFloat(2d + 3)) * miss[d + 2]
        cert = Certificate(equations = N >= 2 ? "P_$N'(x_i) = 0 at the interior Lobatto points (Newton correction)" :
                                                "none (the two end points)",
                           residual = BigFloat(res; precision = 64), residual_bits = wbits,
                           digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)),
                           next_error = BigFloat(next; precision = 64))
        prov = Provenance(family = "Filon", derivation = Derived(),
                          path = ["weight: $(w) on $(dom.base)",
                                  "nodes: the $(N + 1) Gauss–Lobatto points, ±1 and the zeros of P_$N'",
                                  "moments: ∫ P_k e^{iκx} dx = 2iᵏ j_k(κ), κ = ωh = $(Float64(κ)); spherical Bessel " *
                                  "functions by upward recurrence at a precision raised by its loss",
                                  "weights: interpolatory, through the discrete orthogonality of P_k on the Lobatto points",
                                  "degree: $d" * (d > N ? " (interpolatory degree $N, and the kernel's parity)" : " (interpolatory)")],
                          seed_source = "none (derived)", citations = [FILON_1928, ISERLES_NORSETT_2005])
        QuadratureRule(xs, ws, dom, PolynomialDegree(d), prov, cert)
    end
end
