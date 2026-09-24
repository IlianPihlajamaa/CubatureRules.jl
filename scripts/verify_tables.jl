# Verify every shipped symmetric rule, independently of the machinery that produced it
# (PLAN §8: verification as an output, not just a test).
#
#     julia --project -t auto scripts/verify_tables.jl [digits] [degrees]
#
# `degrees` (for example `49:50`) restricts the campaign to those degrees, which is what a
# newly generated seed needs; without it every shipped rule is checked.
#
# For each rule: refine to `digits` (default 60), then check
#   * exactness at its claimed degree and non-exactness one degree higher, against the
#     orthonormal Dubiner basis at twice the precision (`check`)
#   * positive weights, interior nodes, exact symmetry, distinct nodes
#   * every monomial of degree ≤ d integrated to the exact Dirichlet moment, computed in
#     independent arithmetic rather than through the invariant basis
# and, for triangles, compare the point count with the published count of
# Xiao & Gimbutas (2010), Table 1, column n₆.
#
# Exit status is non-zero if any rule fails, so this can be run as a campaign in CI.

using CubatureRules, Printf, TOML
import CubatureRules: degree_range     # public, not exported
const CR = CubatureRules

# Published counts for fully symmetric, positive, interior rules, for comparison only.
# Tetrahedra: Witherden & Vincent (2015) Table 1 (to degree 10), then Zhang, Cui & Liu
# (2009) Table 4.2 (to degree 14).
const TET_PUBLISHED = [1, 4, 8, 14, 14, 24, 35, 46, 59, 81, 109, 140, 171, 236]
const XG_N6 = [1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79,
               87, 96, 103, 112, 120, 130, 141, 150, 159, 171, 181, 193, 204, 214, 228,
               243, 252, 267, 282, 295, 309, 324, 339, 354, 370, 385, 399, 423, 435, 453]

"Largest relative error over every monomial of degree ≤ d, in `bits`-bit arithmetic."
function monomial_error(r, d::Int, bits::Int)
    x = nodes(r)
    w = weights(r)
    D = length(first(x))
    return setprecision(BigFloat, bits) do
        worst = big(0.0)
        for α in CR.compositions(d, D + 1)
            e = α[2:end]                                    # Cartesian exponents
            exact = BigFloat(CR.monomial_moment(Simplex{D}(), Tuple(e)))
            got = sum(BigFloat(w[i]) * prod(BigFloat(x[i][j])^e[j] for j in 1:D) for i in eachindex(x))
            worst = max(worst, abs(got - exact) / exact)
        end
        worst
    end
end

function campaign(digits::Int; degrees = nothing)
    ok = true
    for (fam, dom, published) in ((XiaoGimbutas(), Simplex{2}(), XG_N6), (FullySymmetric(), Simplex{3}(), TET_PUBLISHED))
        name = CR.family_name(fam)
        println("\n", name, " on ", dom, " — refined to ", digits, " digits")
        println("  degree  points  published  exact  sharp  positive  interior  symmetric  min sep   monomials")
        wanted = degrees === nothing ? (1:last(degree_range(fam, dom))) : intersect(degrees, 1:last(degree_range(fam, dom)))
        isempty(wanted) && (println("  (no requested degrees in range)"); continue)
        for d in wanted
            r = rule(fam, dom; degree = d, digits)
            v = check(r)
            x = nodes(r)
            sep = length(x) == 1 ? Inf :
                  minimum(maximum(abs, x[i] - x[j]) for i in eachindex(x) for j in (i + 1):length(x))
            mono = monomial_error(r, d, 4 * CR.digits_to_bits(digits))
            pub = published === nothing || d > length(published) ? "—" : string(published[d])
            good = v.exact && v.sharp !== false && v.positive && v.interior && v.symmetric === true &&
                   sep > 1e-6 && mono < big(10.0)^(-digits + 2)
            ok &= good
            @printf("  %6d  %6d  %9s  %5s  %5s  %8s  %8s  %9s  %.1e  %.1e%s\n", d, npoints(r), pub,
                    v.exact, v.sharp, v.positive, v.interior, v.symmetric, sep, mono, good ? "" : "   <-- FAILED")
            flush(stdout)
        end
    end
    # where the shipped triangle counts differ from the published ones
    println("\nPoint counts differing from Xiao & Gimbutas (2010), Table 1:")
    for d in 1:min(last(degree_range(XiaoGimbutas(), Simplex{2}())), length(XG_N6))
        n = npoints(XiaoGimbutas(), Simplex{2}(), d)
        e = CR.xg_entry_for(d)
        e.degree == d || continue
        n == XG_N6[d] && continue
        @printf("  degree %2d: %3d here, %3d published (%+d)\n", d, n, XG_N6[d], n - XG_N6[d])
    end
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    digits = isempty(ARGS) ? 60 : parse(Int, ARGS[1])
    degrees = length(ARGS) < 2 ? nothing : (:)(parse.(Int, split(ARGS[2], ":"))...)
    exit(campaign(digits; degrees) ? 0 : 1)
end
