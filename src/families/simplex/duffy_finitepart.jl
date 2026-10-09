# Rules for the strongly singular and hypersingular kernels of boundary-element methods on a
# flat triangle, x₀ on it (notes/singular-bem.md, stage 5): the finite parts with respect to
# r = |y − x₀| of
#
#     ∫_T f(y) / r³ dy                  (InverseDistanceCubed)
#     ∫_T f(y) (y − x₀)·e / r³ dy       (InverseDistanceGradient)
#
# On a sub-triangle (x₀, p, q), the Duffy map y = x₀ + s v(t), v = a + t (b − a), gives
# r = s √q(t), q = |v|², and dy = s |a × b| ds dt; with f(x₀ + s v(t)) = Σ c_k(t) sᵏ,
#
#     1/r³:      |a × b| q^{-3/2} f/s² ds dt,        (v·e) kernel:  (v·e) |a × b| q^{-3/2} f/s ds dt.
#
# Excluding r < ε is excluding s < ε/√q(t), so the finite part in r is the finite part in s plus
# what the change of the cut-off leaves: from ∫_{ε/√q}^1 f/s² ds = FP + c₀ √q/ε − c₁ (ln ε − ln √q)
# and ∫_{ε/√q}^1 f/s ds = FP − c₀ (ln ε − ln √q), dropping the terms in ε,
#
#     1/r³:      ∫ dt |a × b| q^{-3/2} [ −c₀ + c₁ ln √q + Σ_{k≥2} c_k/(k − 1) ]
#     (v·e):     ∫ dt (v·e) |a × b| q^{-3/2} [ c₀ ln √q + Σ_{k≥1} c_k/k ]
#
# (Guiggiani's correction of the finite part between local and distance variables). In s, an
# interpolatory rule on d + 1 Gauss–Legendre nodes for the polynomial part of the bracket. The
# log terms are φ(x₀) ∫ (v·e) ln √q q^{-3/2} dt and ∇φ(x₀)·∫ v ln √q q^{-3/2} dt, constants per
# sub-triangle; c₀ = φ(x₀) and c₁(t) = ∇φ(x₀)·v(t) come from the nodes of any ray by
# extrapolation and differentiation weights, so they are folded into the weights and no node
# sits at x₀. In t, the Gauss rule of q^{-3/2}, from the Stieltjes procedure on the
# sinh-substituted discretisation of DuffyGauss (the weight is smooth in the new variable).

"""
    DuffyFinitePart()

Rules for the finite parts of `∫_T f(y) K(y) dy` with the strongly singular kernel
[`InverseDistanceGradient`](@ref) and the hypersingular kernel [`InverseDistanceCubed`](@ref)
on a triangle in the plane or in space, with `x₀` on the triangle: at a vertex, on an edge or
inside. The finite part is taken with respect to the distance `r = |y − x₀|`, so the finite
parts over the panels around a collocation point add up to the one over their union, and for
the gradient kernel with `x₀` inside it is the Cauchy principal value.

The triangle is cut at `x₀`; on each piece, the Duffy map with `d + 1` Gauss–Legendre nodes in
the radial direction, weighted for the finite part, and the Gauss rule of `q(t)^{-3/2}` in the
angular one, with the logarithmic terms that the finite part in the Duffy variable and in `r`
differ by folded into the weights. Exact on polynomials `f` of degree `d`, with
`(d + 1) · ⌈(d + 1)/2⌉` points per piece (at least two rays) for `1/r³` and
`(d + 1) · ⌈(d + 2)/2⌉` for the gradient kernel; the weights are signed, the nodes interior.

```julia
dom = WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistanceCubed((1//4, 1//4)))
r = rule(dom; degree = 6)
integrate(y -> 1 + y[1] * y[2], r)          # ⨎_T (1 + y₁y₂)/|y − x₀|³ dy
```
"""
struct DuffyFinitePart <: RuleFamily end

derivation(::Type{DuffyFinitePart}) = Derived()
family_name(::DuffyFinitePart) = "DuffyFinitePart"

const GUIGGIANI_1992 = Citation(key = "Guiggiani1992",
                                authors = ["M. Guiggiani", "G. Krishnasamy", "T. J. Rudolphi", "F. J. Rizzo"],
                                title = "A general algorithm for the numerical solution of hypersingular boundary integral equations",
                                journal = "Journal of Applied Mechanics", year = 1992, volume = "59", pages = "604–614")

# rays per sub-triangle: the Gauss rule in t must be exact to degree d (d + 1 with the factor
# v·e); 1/r³ needs two rays at least, to recover ∇φ(x₀) in the plane
finitepart_rays(dom, degree) = dom.weight isa InverseDistanceCubed ? max(2, cld(degree + 1, 2)) : max(1, cld(degree + 2, 2))

candidates(::Type{DuffyFinitePart}, dom::FinitePartDomain, ::PolynomialDegree) =
    kernel_geometry(dom).place === :near ? DuffyFinitePart[] : [DuffyFinitePart()]
npoints(::DuffyFinitePart, dom, degree::Integer) = npatches(dom) * (degree + 1) * finitepart_rays(dom, degree)
claimed_degree(::DuffyFinitePart, dom, degree) = Int(degree)
degree_range(::DuffyFinitePart, dom) = 0:typemax(Int)
properties(::DuffyFinitePart, dom, degree) = (positive = false, interior = true, symmetry = :none, nested = false)

"`∫₀¹ g(t) q(t)^{-3/2} dt` discretised by `t = t* + η sinh u` and Gauss–Legendre `(x, w)` in `u`."
function inverse_cube_grid(A, B, (x, w))
    D = B - A
    dd = dot(D, D)
    ts = -dot(A, D) / dd
    η = _area2(A, B) / dd
    u0, u1 = asinh(-ts / η), asinh((1 - ts) / η)
    c, h = (u0 + u1) / 2, (u1 - u0) / 2
    u = [c + h * xi for xi in x]
    # dt = η cosh u du and q = |D|² η² cosh² u, so dt/q^{3/2} = du/(|D|³ η² cosh² u)
    return [ts + η * sinh(ui) for ui in u], [h * wi / (sqrt(dd)^3 * η^2 * cosh(ui)^2) for (ui, wi) in zip(u, w)]
end

"""
    inverse_cube_quadratic_rule(a, b, m, bits, gl) -> (t, λ, M, gap, wbits)

The `m`-point Gauss rule on `[0, 1]` of the weight `|a + t (b − a)|⁻³`, from the Stieltjes
procedure on the sinh-substituted discretisation (as `inverse_sqrt_quadratic_rule`).
"""
function inverse_cube_quadratic_rule(a, b, m::Integer, bits::Integer, gl)
    discretize(M, wbits) = with_bits(wbits) do
        inverse_cube_grid(SVector{length(a),BigFloat}(a), SVector{length(b),BigFloat}(b), gl(M, wbits))
    end
    return converged_stieltjes(discretize, m, bits; name = "DuffyFinitePart", failure = (M, gap) ->
        "the Gauss rule of q(t)^(-3/2) on a sub-triangle with edge vectors $(Float64.(a)), $(Float64.(b)) did not " *
        "converge to $bits bits with $M points (the recurrence still changed by $gap)")
end

"`∫₀¹ g(t) q(t)^{-3/2} dt` to `wbits`, the discretisation doubled until two agree; `g` may be vector-valued."
function inverse_cube_integral(A, B, g, wbits, gl)
    prev = nothing
    M = 64
    while M <= 8192
        t, w = inverse_cube_grid(A, B, gl(M, wbits))
        val = sum(wi * g(ti) for (ti, wi) in zip(t, w))
        prev !== nothing && maximum(abs, val - prev) <= ldexp(one(BigFloat), -(wbits - 24)) * max(1, maximum(abs, val)) &&
            return val
        prev = val
        M *= 2
    end
    throw(RefinementError("DuffyFinitePart", "a logarithmic moment on the sub-triangle with edge vectors " *
                                             "$(Float64.(A)), $(Float64.(B)) did not converge with 8192 points"))
end

"Coefficients `κ` with `Σ κⱼ vⱼ = L` for vectors `vⱼ` and `L` in the plane of `a` and `b` (minimum norm)."
function inplane_coefficients(a, b, vs, L)
    if length(a) == 2
        Vm = reduce(hcat, vs)
        return transpose(Vm) * ((Vm * transpose(Vm)) \ L)
    end
    e1 = a / norm(a)
    e2 = b - dot(b, e1) * e1
    e2 /= norm(e2)
    Vm = reduce(hcat, [[dot(v, e1), dot(v, e2)] for v in vs])
    return transpose(Vm) * ((Vm * transpose(Vm)) \ [dot(L, e1), dot(L, e2)])
end

function build(f::DuffyFinitePart, dom::FinitePartDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("nodes for the finite-part kernels are irrational; $(T) is not supported"))
    g = kernel_geometry(dom)
    g.place === :near && throw(ArgumentError("DuffyFinitePart is for x₀ on the triangle; $(Tuple(dom.weight.x0)) is off it"))
    x0, patches, place, moved = g.x0, g.patches, g.place, g.moved
    D = length(x0)
    hyper = dom.weight isa InverseDistanceCubed
    ms, mt = degree + 1, finitepart_rays(dom, degree)
    # the interpolatory weights in s grow with the degree; the logarithmic moments are
    # integrals of smooth functions: a guard for both
    wbits = ctx.bits + 8degree + 48
    cache = Dict{Tuple{Int,Int},Tuple{Vector{BigFloat},Vector{BigFloat}}}()
    gl(M, bits) = get!(() -> gauss_jacobi_work(M, 0, 0, bits)[1:2], cache, (M, bits))
    trules = [inverse_cube_quadratic_rule(p - x0, q - x0, mt, wbits, gl) for (p, q) in patches]
    xs, ws = SVector{D,T}[], T[]
    with_bits(wbits) do
        sx, _ = gl(ms, wbits)
        s = (1 .+ sx) ./ 2
        Vs = [s[i]^k for k in 0:(ms - 1), i in 1:ms]
        weights_for(μ) = Vs \ [μ(k) for k in 0:(ms - 1)]
        # the finite part in s of the polynomial part, and f(0) or f'(0), from the nodes of a ray
        ω = hyper ? weights_for(k -> k == 0 ? -one(BigFloat) : k == 1 ? zero(BigFloat) : 1 / BigFloat(k - 1)) :
                    weights_for(k -> k == 0 ? zero(BigFloat) : 1 / BigFloat(k))
        χ = weights_for(k -> k == (hyper ? 1 : 0) ? one(BigFloat) : zero(BigFloat))
        X0 = SVector{D,BigFloat}(x0)
        e = hyper ? nothing : SVector{D,BigFloat}(map(BigFloat, dom.weight.e))
        for ((p, q), (t, λ, _, _, _)) in zip(patches, trules)
            A, B = SVector{D,BigFloat}(p - x0), SVector{D,BigFloat}(q - x0)
            J = _area2(A, B)
            V(tj) = A + tj * (B - A)
            # the logarithmic term, as coefficients on the rays: Σ_j κ_j c_1(t_j) = ∇φ(x₀)·L for
            # 1/r³, Σ_j κ_j c_0 = φ(x₀) C for the gradient kernel
            κ = if hyper
                inplane_coefficients(A, B, V.(t), inverse_cube_integral(A, B, tj -> V(tj) * log(norm(V(tj))), wbits, gl))
            else
                inverse_cube_integral(A, B, tj -> dot(V(tj), e) * log(norm(V(tj))), wbits, gl) .* λ ./ sum(λ)
            end
            for (j, tj) in enumerate(t), i in eachindex(s)
                y = X0 + s[i] * V(tj)
                radial = hyper ? λ[j] : dot(V(tj), e) * λ[j]
                push!(xs, SVector{D,T}(ntuple(k -> finalize_number(ctx, y[k]), D)))
                push!(ws, finalize_number(ctx, J * (radial * ω[i] + κ[j] * χ[i])))
            end
        end
    end
    gap = maximum(r[4] for r in trules)
    kname = hyper ? "1/|y − x₀|³" : "(y − x₀)·e/|y − x₀|³"
    cert = Certificate(equations = "recurrence of the discretised weight q(t)^(-3/2) on each sub-triangle: " *
                                   "$(join((r[3] for r in trules), ", ")) against half as many points",
                       residual = BigFloat(gap; precision = 64), residual_bits = wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
    prov = Provenance(family = "DuffyFinitePart", derivation = Derived(),
                      path = ["kernel: finite part (with respect to r = |y − x₀|) of $kname on $(dom.base), x₀ = " *
                              "$(Tuple(dom.weight.x0)) " *
                              (place === :vertex ? "a vertex" : place === :edge ? "on an edge" : "inside") *
                              ": $(length(patches)) sub-triangle$(length(patches) == 1 ? "" : "s") with a vertex at x₀" *
                              (moved > 0 ? @sprintf("; x₀ moved by %.1e onto the triangle to absorb rounding, so the rule is for x₀ = %s",
                                                    moved, Tuple(Float64.(x0))) : ""),
                              "Duffy map y = x₀ + s (a + t (b − a)) on each, r = s √q(t)",
                              "s: $ms Gauss–Legendre nodes on [0, 1], interpolatory weights for the finite part in s",
                              "t: $mt-point Gauss rule of q(t)^(-3/2), from the Stieltjes procedure on a sinh-substituted " *
                              "discretisation",
                              "the terms by which the finite parts in s and in r differ, " *
                              (hyper ? "∇φ(x₀)·∫ v ln √q q^(-3/2) dt" : "φ(x₀) ∫ (v·e) ln √q q^(-3/2) dt") *
                              ", folded into the weights by " * (hyper ? "differentiation" : "extrapolation") *
                              " to s = 0 along the rays"],
                      seed_source = "none (derived)", citations = [DUFFY_1982, GUIGGIANI_1992, GAUTSCHI_1982])
    return QuadratureRule(xs, ws, dom, PolynomialDegree(degree), prov, cert)
end
