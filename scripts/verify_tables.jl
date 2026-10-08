# Verify every shipped symmetric rule, independently of the machinery that produced it
# (PLAN §8: verification as an output, not just a test).
#
#     julia --project -t auto scripts/verify_tables.jl [digits] [degrees] [domains]
#
# `degrees` (for example `49:50`, or `11,22`) restricts the campaign to those degrees, and
# `domains` (any of `triangle,tetrahedron,simplex4,square,cube,hypercube,disk,ball,pyramid,wedge`)
# to those domains, which is what newly generated seeds need; without them every shipped
# rule is checked.
#
# For each rule: refine to `digits` (default 60), then check
#   * exactness at its claimed degree and non-exactness one degree higher, against the
#     orthonormal Dubiner basis at twice the precision (`check`)
#   * positive weights, interior nodes, exact symmetry, distinct nodes
#   * every monomial of degree ≤ d integrated to its exact moment (Dirichlet on the simplex,
#     products of 2/(k+1) on the box, Gamma functions on the ball, beta integrals on the
#     pyramid and the wedge), computed in independent arithmetic rather than through the
#     invariant basis or the orthonormal one
# and compare the point count with the published one where there is one: Xiao & Gimbutas
# (2010), Table 1, column n₆, for triangles, and the tables below for tetrahedra, boxes,
# pyramids and wedges. Simplices, pyramids and wedges are checked at every degree; boxes,
# disks and balls at the odd degrees, which are the ones tabulated.
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

# Witherden & Vincent (2015): fully symmetric quadrilateral and hexahedron rules
const SQUARE_PUBLISHED = Dict(1 => 1, 3 => 4, 5 => 8, 7 => 12, 9 => 20, 11 => 28, 13 => 37, 15 => 48,
                              17 => 60, 19 => 72, 21 => 85)
const CUBE_PUBLISHED = Dict(1 => 1, 3 => 6, 5 => 14, 7 => 34, 9 => 58, 11 => 90)
const PYRAMID_PUBLISHED = Dict(1 => 1, 2 => 5, 3 => 6, 4 => 10, 5 => 15, 6 => 24, 7 => 31, 8 => 47, 9 => 62,
                               10 => 83)
const WEDGE_PUBLISHED = Dict(1 => 1, 2 => 5, 3 => 8, 4 => 11, 5 => 16, 6 => 28, 7 => 35, 8 => 46, 9 => 60, 10 => 85)

# ∫ |x^e|: the moment itself on the simplex, where every one is positive; on the box, where
# the odd ones vanish, the moment of |x^e|, so that the error is relative to the integrand
moment_scale(::Simplex{D}, e) where {D} = CR.monomial_moment(Simplex{D}(), Tuple(e))
moment_scale(::Orthotope, e) = prod(big(2) // (k + 1) for k in e)
# on the ball |xᵉ| ≥ |x|^(e + 1) for the odd exponents, so the moment of the next even
# exponent is a lower bound for ∫ |xᵉ|, and the error is measured against at most that
moment_scale(::Ball{D}, e) where {D} = CR.monomial_moment(Ball{D}(), Tuple(k + isodd(k) for k in e))
# on the pyramid ∫ |xᵉ| is the moment with the parity of the base exponents ignored
moment_scale(::Pyramid, e) = big(4) // ((e[1] + 1) * (e[2] + 1)) * factorial(big(e[3])) *
                             factorial(big(e[1] + e[2] + 2)) // factorial(big(sum(e) + 3))
# on the wedge the triangle's moment times ∫ |z|^c
moment_scale(::Wedge, e) = CR.barycentric_moment((0, e[1], e[2])) * big(2) // (e[3] + 1)

"Largest relative error over every monomial of degree ≤ d, in `bits`-bit arithmetic."
function monomial_error(r, dom, d::Int, bits::Int)
    x = nodes(r)
    w = weights(r)
    D = length(first(x))
    return setprecision(BigFloat, bits) do
        worst = big(0.0)
        for α in CR.compositions(d, D + 1)
            e = α[2:end]                                    # Cartesian exponents
            exact = BigFloat(CR.monomial_moment(dom, Tuple(e)))
            got = sum(BigFloat(w[i]) * prod(BigFloat(x[i][j])^e[j] for j in 1:D) for i in eachindex(x))
            worst = max(worst, abs(got - exact) / BigFloat(moment_scale(dom, e)))
        end
        worst
    end
end

function campaign(digits::Int; degrees = nothing, domains = nothing)
    ok = true
    for (label, fam, dom, published) in (("triangle", XiaoGimbutas(), Simplex{2}(), XG_N6),
                                         ("tetrahedron", FullySymmetric(), Simplex{3}(), TET_PUBLISHED),
                                         ("simplex4", FullySymmetric(), Simplex{4}(), Dict{Int,Int}()),
                                         ("square", FullySymmetric(), Orthotope{2}(), SQUARE_PUBLISHED),
                                         ("cube", FullySymmetric(), Orthotope{3}(), CUBE_PUBLISHED),
                                         ("hypercube", FullySymmetric(), Orthotope{4}(), Dict{Int,Int}()),
                                         ("disk", FullySymmetric(), Disk(), Dict{Int,Int}()),
                                         ("ball", FullySymmetric(), Ball{3}(), Dict{Int,Int}()),
                                         ("pyramid", FullySymmetric(), Pyramid(), PYRAMID_PUBLISHED),
                                         ("wedge", FullySymmetric(), Wedge(), WEDGE_PUBLISHED))
        domains === nothing || label in domains || continue
        name = CR.family_name(fam)
        println("\n", name, " on ", dom, " — refined to ", digits, " digits")
        println("  degree  points  published  exact  sharp  positive  interior  symmetric  min sep   monomials")
        top = last(degree_range(fam, dom))
        shipped = dom isa Union{Orthotope,Ball} ? (1:2:top) : (1:top)
        wanted = degrees === nothing ? shipped : intersect(degrees, shipped)
        isempty(wanted) && (println("  (no requested degrees in range)"); continue)
        for d in wanted
            r = rule(fam, dom; degree = d, digits)
            v = check(r)
            x = nodes(r)
            sep = length(x) == 1 ? Inf :
                  minimum(maximum(abs, x[i] - x[j]) for i in eachindex(x) for j in (i + 1):length(x))
            mono = monomial_error(r, dom, d, 4 * CR.digits_to_bits(digits))
            p = get(published, d, nothing)
            pub = p === nothing ? "—" : string(p)
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
    degrees = length(ARGS) < 2 ? nothing :
              occursin(":", ARGS[2]) ? (:)(parse.(Int, split(ARGS[2], ":"))...) : parse.(Int, split(ARGS[2], ","))
    domains = length(ARGS) < 3 ? nothing : split(ARGS[3], ",")
    exit(campaign(digits; degrees, domains) ? 0 : 1)
end
