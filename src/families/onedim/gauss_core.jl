# The Gauss driver shared by every classical family (PLAN §1, §6 Tier 1).
#
# A Gauss rule is fixed by the three-term recurrence of the orthonormal polynomials for its
# weight:
#
#     b_{k+1} p_{k+1}(x) = (x - a_k) p_k(x) - b_k p_{k-1}(x),    μ₀ = ∫ w
#
# so a family only has to supply a_k, b_k and μ₀. The pipeline is the same in each case:
# Golub–Welsch in Float64 for the seed, Newton on the recurrence at working precision, and
# weights from the Christoffel function w_i = μ₀ / Σ_k q_k(x_i)².

"Guard bits for Newton on a three-term recurrence: O(log n), which is right there (§2.7)."
gj_guard_bits(n) = 24 + 2 * ceil(Int, log2(n + 1))

"""
    Recurrence(a, b, mass)

The orthonormal three-term recurrence of a weight: `a(k, T)` and `b(k, T)` give the
coefficients in type `T`, and `mass(T)` gives `μ₀ = ∫ w`.
"""
struct Recurrence{A,B,M}
    a::A
    b::B
    mass::M
end

"""
    gauss_from_recurrence(n, rec, bits) -> (x, w, iterations)

`n`-point Gauss nodes (ascending) and weights for the weight described by `rec`, as
BigFloats of precision `bits`. Callers choose `bits` including any guard; nothing here
reads the ambient precision.
"""
function gauss_from_recurrence(n::Integer, rec::Recurrence, bits::Integer; name::String = "Gauss")
    n >= 1 || throw(ArgumentError("need at least one node"))
    # --- seed: Golub–Welsch in Float64
    a64 = [Float64(rec.a(k, Float64)) for k in 0:(n - 1)]
    b64 = [Float64(rec.b(k, Float64)) for k in 1:(n - 1)]
    seed = n == 1 ? a64 : eigvals(SymTridiagonal(a64, b64))
    sort!(seed)
    return with_bits(bits) do
        a = [rec.a(k, BigFloat) for k in 0:n]
        b = [rec.b(k, BigFloat) for k in 1:n]
        μ = BigFloat(rec.mass(BigFloat))
        x = BigFloat.(seed)
        tol = ldexp(BigFloat(1), -(bits - 6))
        # The rounding noise in q/q' grows with the size of the recurrence coefficients —
        # roughly the largest node — not with the node being refined. At Laguerre n = 300 it
        # sat at 160 ε on a node near 0.025, above a tolerance of 64 ε relative to max(|x|, 1),
        # and Newton, converged, never met its stopping test. The guard bits cover the looser
        # absolute tolerance this gives the smallest nodes.
        scale = max(one(BigFloat), BigFloat(maximum(abs, seed)))
        iters = 0
        for i in 1:n
            xi = x[i]
            last_small = false
            for it in 1:100
                q, dq, _ = gauss_eval(n, xi, a, b)
                δ = q / dq
                xi -= δ
                iters = max(iters, it)
                small = abs(δ) <= tol * max(abs(xi), scale)
                small && last_small && break
                last_small = small
                it == 100 && throw(RefinementError(name, "Newton on p_$n did not converge at node $i"))
            end
            x[i] = xi
        end
        w = similar(x)
        for i in 1:n
            _, _, s = gauss_eval(n, x[i], a, b)
            w[i] = μ / s
        end
        x, w, iters
    end
end

# Unnormalised (q₀ = 1) orthonormal recurrence: returns q_n(x), q_n'(x) and Σ_{k<n} q_k(x)².
function gauss_eval(n, x, a, b)
    q0, q1 = zero(x), one(x)
    d0, d1 = zero(x), zero(x)
    s = one(x)
    for k in 0:(n - 1)
        bk = k == 0 ? zero(x) : b[k]
        q2 = ((x - a[k + 1]) * q1 - bk * q0) / b[k + 1]
        d2 = (q1 + (x - a[k + 1]) * d1 - bk * d0) / b[k + 1]
        q0, q1 = q1, q2
        d0, d1 = d1, d2
        k < n - 1 && (s += q1^2)
    end
    return q1, d1, s
end

"Newton-correction residual `max |q_n(x̂)/q_n'(x̂)|` of the rounded nodes, evaluated at `bits`."
function gauss_residual(xhat, rec::Recurrence, bits)
    n = length(xhat)
    with_bits(bits) do
        a = [rec.a(k, BigFloat) for k in 0:n]
        b = [rec.b(k, BigFloat) for k in 1:n]
        maximum(xhat) do xi
            q, dq, _ = gauss_eval(n, BigFloat(xi), a, b)
            abs(q / dq)
        end
    end
end
