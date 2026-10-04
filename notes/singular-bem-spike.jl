# Feasibility spike for vertex-singular triangle rules, ∫_T p(y) / |y − v₀| dy, for
# notes/singular-bem.md. Run with `julia --project notes/singular-bem-spike.jl` (about a minute).
#
# Duffy from v₀: y = v₀ + s (e₁ + t (e₂ − e₁)), dy = s |e₁ × e₂| ds dt, |y − v₀| = s √q(t),
# q(t) = |e₁ + t (e₂ − e₁)|². So p(y)/|y − v₀| dy = p(s, t) |e₁ × e₂| / √q(t) ds dt, with p of
# degree ≤ d in s and in t: Gauss–Legendre in s times the Gauss rule of the weight 1/√q(t)
# should be exact. The reference is independent: polar coordinates about v₀, the radial
# integral in closed form, the angular one by Gauss–Legendre at high precision.
using CubatureRules, Printf, LinearAlgebra
const CR = CubatureRules

const DIGITS = 40
setprecision(BigFloat, 4 * CR.digits_to_bits(DIGITS))

v0 = big.([0.0, 0.0]); v1 = big.([1.0, 0.0]); v2 = [big"0.3", big"0.8"]
e1, e2 = v1 - v0, v2 - v0
cross2(a, b) = a[1] * b[2] - a[2] * b[1]
area2 = abs(cross2(e1, e2))
q(t) = sum(abs2, e1 + t * (e2 - e1))

# reference: ∫_T y^α / |y| dy = ∫_{θ₀}^{θ₁} cos^a θ sin^b θ R(θ)^{|α|+1} / (|α|+1) dθ, with R(θ)
# the distance from v₀ to the opposite edge along direction θ (v₀ is the origin here)
θ0, θ1 = atan(e1[2], e1[1]), atan(e2[2], e2[1])
R(θ) = cross2(v1, v2 - v1) / cross2([cos(θ), sin(θ)], v2 - v1)
gl = rule(Interval(); degree = 399, digits = 3 * DIGITS)          # 200 points, mapped below
function reference(a, b)
    k = a + b
    sum(w * (θ1 - θ0) / 2 * (cos(θ)^a * sin(θ)^b * R(θ)^(k + 1) / (k + 1))
        for (x, w) in zip(nodes(gl), weights(gl)) for θ in ((θ0 + θ1) / 2 + (θ1 - θ0) / 2 * x,))
end

function duffy_rule(d; weighted = true)
    rs = rule(Interval(0, 1); degree = d, digits = DIGITS)
    rt = weighted ? rule(WeightedDomain(Interval(0, 1), FunctionWeight(t -> 1 / sqrt(q(t)), 0, 1)); degree = d, digits = DIGITS) :
                    rule(Interval(0, 1); degree = d, digits = DIGITS)
    pts, wts = Vector{BigFloat}[], BigFloat[]
    for (s, ws) in zip(nodes(rs), weights(rs)), (t, wt) in zip(nodes(rt), weights(rt))
        y = v0 + s * (e1 + t * (e2 - e1))
        # weighted: the 1/√q(t) is in the t-rule; unweighted: the kernel is in the integrand
        push!(pts, y); push!(wts, ws * wt * area2 * (weighted ? 1 : 1 / sqrt(q(t))))
    end
    return pts, wts
end

errs(pts, wts, d) = maximum(abs(sum(w * y[1]^a * y[2]^b for (y, w) in zip(pts, wts)) - reference(a, b)) / abs(reference(a, b))
                             for a in 0:d for b in 0:(d - a))

println("vertex-singular rule, max relative error over y^α/|y| with |α| ≤ d (and at d + 1):")
for d in (2, 5, 10, 20)
    t = @elapsed pts, wts = duffy_rule(d)
    e, e1_ = errs(pts, wts, d), errs(pts, wts, d + 1)
    pu, wu = duffy_rule(d; weighted = false)
    @printf("  d = %2d: %3d points, error %.1e (at d+1: %.1e)  [%.1f s]; Duffy × Gauss–Legendre: %.1e\n",
            d, length(wts), Float64(e), Float64(e1_), t, Float64(errs(pu, wu, d)))
    flush(stdout)
end

# a standard (non-singular) triangle rule of the same degree, for comparison: the kernel is
# not polynomial, so it converges only algebraically
T = Simplex([0, 0], [1, 0], [3//10, 8//10])
for d in (10, 20, 40)
    r = map_to(rule(Simplex{2}(); degree = d, digits = DIGITS), T)
    got = sum(w / sqrt(sum(abs2, x)) for (x, w) in zip(nodes(r), weights(r)))
    @printf("  smooth rule of degree %2d (%3d points) on 1/|y|: relative error %.1e\n", d, npoints(r),
            Float64(abs(got - reference(0, 0)) / reference(0, 0)))
end

# Helmholtz e^{ikr}/r: the rule is exact for p/r, and e^{ikr}/r = Σ (ik)^n r^{n-1}/n! has the
# smooth terms r^{2m} too, on which the 1/√q rule is not exact — how far does it get?
k = big(5)
refH = let
    # ∫ e^{ikr}/r dy in polar coordinates: ∫ dθ ∫_0^R e^{ikr} dr = ∫ (e^{ikR} − 1)/(ik) dθ
    sum(w * (θ1 - θ0) / 2 * (exp(im * k * R(θ)) - 1) / (im * k)
        for (x, w) in zip(nodes(gl), weights(gl)) for θ in ((θ0 + θ1) / 2 + (θ1 - θ0) / 2 * x,))
end
println("Helmholtz e^{ikr}/r, k = 5:")
for d in (5, 10, 20, 30)
    pts, wts = duffy_rule(d)
    # the rule's weights already carry 1/|y|: integrate e^{ikr} against them
    got = sum(w * exp(im * k * norm(y)) for (y, w) in zip(pts, wts))
    @printf("  d = %2d: %3d points, relative error %.1e\n", d, length(wts), Float64(abs(got - refH) / abs(refH)))
end
