# Rules for the kernel 1/|y − x₀| on a tetrahedron with x₀ in it: the weakly singular volume
# integrals of volume integral equations (notes/singular-bem.md, stage 6).
#
# The tetrahedron is the union of the cones from x₀ over its faces (those whose planes do not
# contain x₀). On the cone over a face F at the distance h from x₀, y = x₀ + s (z − x₀) with z
# on F and s in [0, 1] gives dy = h s² ds dA_z and |y − x₀| = s |z − x₀|, so
#
#     f(y)/|y − x₀| dy = h · s f(x₀ + s (z − x₀)) ds · dA_z/|z − x₀|.
#
# For a polynomial f of degree d that is a polynomial of degree ≤ d + 1 in s, which
# Gauss–Legendre integrates exactly, and for each s a polynomial of degree ≤ d in z against
# 1/|z − x₀| on F with x₀ off F's plane: the near-singular kernel on a triangle in space,
# which DuffySinh integrates exactly to working precision however close x₀ is to the face.

"""
    DuffyCone()

Rules for `∫_K f(y)/|y − x₀| dy` on a tetrahedron `K` with `x₀` in it (inside, on a face, on
an edge or at a vertex), the weakly singular kernel of volume integral equations: an
[`InverseDistance`](@ref) weight on a `Simplex{3}` with three-dimensional vertices.

The tetrahedron is cut into the cones from `x₀` over its faces; on the cone over a face at the
distance `h`, `y = x₀ + s (z − x₀)` turns `f(y)/|y − x₀| dy` into `h s f ds` times
`dA_z/|z − x₀|` on the face. Gauss–Legendre in `s` is exact, and on the face, where `x₀` is off
the plane, the [`DuffySinh`](@ref) rule for that near-singular kernel. Positive weights,
interior nodes; like `DuffySinh`, the size grows with the precision asked for, and `npoints`
gives the smallest.

```julia
K = Simplex((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
r = rule(WeightedDomain(K, InverseDistance((1//4, 1//5, 1//6))); degree = 5)
integrate(y -> 1 + y[3], r)                 # ∫_K (1 + y₃)/|y − x₀| dy
```
"""
struct DuffyCone <: RuleFamily end

derivation(::Type{DuffyCone}) = Derived()
family_name(::DuffyCone) = "DuffyCone"

candidates(::Type{DuffyCone}, dom::InverseDistanceTetDomain, ::PolynomialDegree) =
    tet_kernel_geometry(dom).inside ? [DuffyCone()] : DuffyCone[]
npoints(::DuffyCone, dom, degree::Integer) =
    duffy_m(degree + 1) * sum((npoints(DuffySinh(), fd, degree) for (fd, _) in cone_faces(dom, Float64)); init = 0)
claimed_degree(::DuffyCone, dom, degree) = min(2duffy_m(degree) - 1, 2duffy_m(degree + 1) - 2)
degree_range(::DuffyCone, dom) = 0:typemax(Int)
properties(::DuffyCone, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

function build(f::DuffyCone, dom::InverseDistanceTetDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("nodes for the kernel 1/|y − x₀| are irrational; $(T) is not supported"))
    g = tet_kernel_geometry(dom)
    g.inside || throw(ArgumentError("DuffyCone is for x₀ in the tetrahedron; $(Tuple(dom.weight.x0)) is outside it"))
    ms = duffy_m(degree + 1)                       # Gauss–Legendre in s, exact to degree 2ms − 1 ≥ d + 1
    wbits = ctx.bits + 16
    # the face rules to the target precision, no more: the radial direction is exact and the
    # weights positive, so combining them loses nothing, and DuffySinh's size grows with the
    # precision asked of it (BigFloat nodes, so that the combination rounds once)
    fctx = BuildContext{BigFloat}(ctx.bits; cancel = ctx.cancel, verbose = max(ctx.verbose - 1, 0))
    faces = with_bits(() -> cone_faces(dom, BigFloat), wbits)
    frules = [build(DuffySinh(), fd, degree, fctx) for (fd, _) in faces]
    K = min(minimum(r -> exactness(r).d, frules), 2ms - 2)
    ctx.verbose >= 1 && @info @sprintf("DuffyCone: %d cone%s, %d points in s, %s on the faces", length(faces),
                                       length(faces) == 1 ? "" : "s", ms, join(npoints.(frules), ", "))
    xs, ws = SVector{3,T}[], T[]
    with_bits(wbits) do
        sx, sw = gauss_jacobi_work(ms, 0, 0, wbits)[1:2]
        X0 = SVector{3,BigFloat}(g.x0)
        for ((_, h), r) in zip(faces, frules), (z, λ) in zip(nodes(r), weights(r)), (si, σi) in zip(sx, sw)
            s, σ = (1 + si) / 2, σi / 2
            y = X0 + s * (SVector{3,BigFloat}(z) - X0)
            push!(xs, SVector{3,T}(ntuple(k -> finalize_number(ctx, y[k]), 3)))
            push!(ws, finalize_number(ctx, h * σ * s * λ))
        end
    end
    cert = Certificate(equations = "the face rules' check, on each cone: " * frules[1].certificate.equations,
                       residual = maximum(r -> r.certificate.residual, frules), residual_bits = wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
    prov = Provenance(family = "DuffyCone", derivation = Derived(),
                      path = vcat(["kernel: 1/|y − x₀| on $(dom.base), x₀ = $(Tuple(dom.weight.x0)): $(length(faces)) " *
                                   "cone$(length(faces) == 1 ? "" : "s") from x₀ over the faces whose planes miss it" *
                                   (g.moved > 0 ? @sprintf("; x₀ moved by %.1e onto a face to absorb rounding", g.moved) : ""),
                                   "y = x₀ + s (z − x₀): f/|y − x₀| dy = h s f ds dA_z/|z − x₀|",
                                   "s: $ms-point Gauss–Legendre on [0, 1]"],
                                  ["face $i (h = $(@sprintf("%.3g", Float64(h)))): DuffySinh, $(npoints(r)) points"
                                   for (i, ((_, h), r)) in enumerate(zip(faces, frules))]),
                      seed_source = "none (derived)", citations = [DUFFY_1982, JOHNSTON_ELLIOTT_2005])
    return QuadratureRule(xs, ws, dom, PolynomialDegree(K), prov, cert)
end
