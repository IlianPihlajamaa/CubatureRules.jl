# The weakly singular kernel of boundary-element methods on a triangle (PLAN §0.3, §6 Tier 5;
# notes/singular-bem.md):
#
#     ∫_T f(y) / |y − x₀| dy,    x₀ in the closed triangle T,
#
# on a triangle in the plane (`Simplex{2}`) or in space (`SurfaceTriangle`). The kernel is the
# weight of a `WeightedDomain`, so a rule for it is an ordinary rule whose weights carry
# 1/|y − x₀|, exact on polynomials f up to its degree.

"""
    PointKernel{D}

A kernel on a triangle in `D` dimensions that is singular at one point `x0`:
[`InverseDistance`](@ref), [`InverseDistanceCubed`](@ref) and
[`InverseDistanceGradient`](@ref).
"""
abstract type PointKernel{D} end

"""
    InverseDistance(x₀)

The weight `1 / |y − x₀|` on a triangle, the weakly singular kernel of boundary-element
methods, as the weight of a [`WeightedDomain`](@ref) on a triangle in the plane (a `Simplex`
with two-dimensional vertices) or in space (a [`SurfaceTriangle`](@ref)). `x₀` may lie on the
triangle (at a vertex, on an edge or inside: a singular integral) or off it (a near-singular
one, hard for ordinary rules when `x₀` is close).

```julia
T = Simplex((0, 0), (1, 0), (0, 1))
dom = WeightedDomain(T, InverseDistance((1//3, 1//3)))
r = rule(dom; degree = 10)
integrate(y -> 1 + y[1]^2, r)            # ∫_T (1 + y₁²) / |y − x₀| dy
```

Rules come from [`DuffyGauss`](@ref) for `x₀` on the triangle and [`DuffySinh`](@ref) off
it. Their weights carry the kernel, so `integrate(f, r)` needs only the smooth factor `f`, and
the rule is exact when `f` is a polynomial of degree up to the rule's degree.

Exact coordinates (integers and rationals) are taken as given. Floating-point ones carry
rounding, and a point meant to lie on the triangle rarely does exactly: the centroid of a
triangle in space, computed in `Float64`, is off its plane by about an ulp. So `x₀` within 64
units in the last place (of the largest coordinate) of the triangle's plane, or of its
boundary, is moved onto it, and the rule is built for that point; the provenance says how
far it moved. A point further off is taken as given, as a near-singular point.
"""
struct InverseDistance{D,P<:Real} <: PointKernel{D}
    x0::SVector{D,P}
    function InverseDistance{D,P}(x0::SVector{D,P}) where {D,P<:Real}
        all(isfinite, x0) || throw(ArgumentError("the singular point must be finite, got $(Tuple(x0))"))
        return new{D,P}(x0)
    end
end
InverseDistance(x0::SVector{D,P}) where {D,P<:Real} = InverseDistance{D,P}(x0)
InverseDistance(x0::Union{Tuple,AbstractVector}) = InverseDistance(SVector{length(x0)}(promote(x0...)))

Base.show(io::IO, w::InverseDistance) = print(io, "InverseDistance(", Tuple(w.x0), ")")

"""
    InverseDistanceCubed(x₀)

The hypersingular kernel `1 / |y − x₀|³` on a triangle, as the weight of a
[`WeightedDomain`](@ref) on a triangle in the plane or in space with `x₀` on it. The integral
does not exist; what the rule computes is its Hadamard finite part with respect to the
distance `r = |y − x₀|`: the integral over the triangle outside the disk `r < ε`, less its
terms in `1/ε` and `ln ε`, as `ε → 0`. That is additive over triangles, so the finite parts
over the triangles around a collocation point add up to the finite part over their union.

On a flat panel it is the Laplace hypersingular kernel, `∂²/∂n_x∂n_y (1/r) = (n_x·n_y)/r³`,
also at an edge or vertex shared with a panel at an angle (the term in `n_y·(y − x₀)` is zero
on the panel). For Helmholtz, `∂²/∂n_x∂n_y (e^{ikr}/r) = e^{ikr}(1 − ikr)/r³` on a flat panel:
integrate `f(y) e^{ikr}(1 − ikr)`, which is smooth along every ray from `x₀` with no linear
term in `r` to spoil the finite part.

```julia
dom = WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistanceCubed((1//4, 1//4)))
r = rule(dom; degree = 6)
integrate(y -> 1 + y[1], r)                 # ⨎_T (1 + y₁)/|y − x₀|³ dy
```

Rules come from [`DuffyFinitePart`](@ref); their weights are signed. `x₀` is read as for
[`InverseDistance`](@ref).
"""
struct InverseDistanceCubed{D,P<:Real} <: PointKernel{D}
    x0::SVector{D,P}
    function InverseDistanceCubed{D,P}(x0::SVector{D,P}) where {D,P<:Real}
        all(isfinite, x0) || throw(ArgumentError("the singular point must be finite, got $(Tuple(x0))"))
        return new{D,P}(x0)
    end
end
InverseDistanceCubed(x0::SVector{D,P}) where {D,P<:Real} = InverseDistanceCubed{D,P}(x0)
InverseDistanceCubed(x0::Union{Tuple,AbstractVector}) = InverseDistanceCubed(SVector{length(x0)}(promote(x0...)))
Base.show(io::IO, w::InverseDistanceCubed) = print(io, "InverseDistanceCubed(", Tuple(w.x0), ")")

"""
    InverseDistanceGradient(x₀, e)

The strongly singular kernel `(y − x₀)·e / |y − x₀|³ = e·∇_{x₀}(1/|y − x₀|)` on a triangle,
as the weight of a [`WeightedDomain`](@ref) on a triangle in the plane or in space with `x₀`
on it. The rule computes its finite part with respect to `r = |y − x₀|`, as for
[`InverseDistanceCubed`](@ref) (here only a `ln ε` term is dropped); for `x₀` inside the
triangle that is the Cauchy principal value.

With `e` the normal at the collocation point it is the Laplace adjoint double-layer kernel,
`n_x·∇_x (1/r)`: zero on the panel of `x₀`, nonzero on a panel meeting it at an edge or vertex
at an angle (only the component of `e` in the triangle's plane matters). For Helmholtz,
`n_x·∇_x (e^{ikr}/r)` is this kernel times `e^{ikr}(1 − ikr)`: integrate `f(y) e^{ikr}(1 − ikr)`.

Rules come from [`DuffyFinitePart`](@ref); their weights are signed.
"""
struct InverseDistanceGradient{D,P<:Real} <: PointKernel{D}
    x0::SVector{D,P}
    e::SVector{D,P}
    function InverseDistanceGradient{D,P}(x0::SVector{D,P}, e::SVector{D,P}) where {D,P<:Real}
        all(isfinite, x0) || throw(ArgumentError("the singular point must be finite, got $(Tuple(x0))"))
        all(isfinite, e) || throw(ArgumentError("the direction must be finite, got $(Tuple(e))"))
        return new{D,P}(x0, e)
    end
end
function InverseDistanceGradient(x0::Union{Tuple,AbstractVector,SVector}, e::Union{Tuple,AbstractVector,SVector})
    length(x0) == length(e) || throw(DimensionMismatch("x₀ and e need the same dimension"))
    v = promote(x0..., e...)
    D = length(x0)
    return InverseDistanceGradient{D,eltype(v)}(SVector{D}(v[1:D]), SVector{D}(v[(D + 1):end]))
end
Base.show(io::IO, w::InverseDistanceGradient) = print(io, "InverseDistanceGradient(", Tuple(w.x0), ", ", Tuple(w.e), ")")

"""
    InverseDistanceDomain

A triangle in the plane (`Simplex{2}`) or in space ([`SurfaceTriangle`](@ref)) carrying an
[`InverseDistance`](@ref) kernel.
"""
const InverseDistanceDomain = Union{WeightedDomain{2,<:Any,<:Simplex{2},<:InverseDistance{2}},
                                    WeightedDomain{3,<:Any,<:SurfaceTriangle,<:InverseDistance{3}}}

"A triangle in the plane or in space carrying a [`PointKernel`](@ref)."
const PointKernelDomain = Union{WeightedDomain{2,<:Any,<:Simplex{2},<:PointKernel{2}},
                                WeightedDomain{3,<:Any,<:SurfaceTriangle,<:PointKernel{3}}}

"A triangle carrying one of the finite-part kernels, [`InverseDistanceCubed`](@ref) or [`InverseDistanceGradient`](@ref)."
const FinitePartDomain = Union{
    WeightedDomain{2,<:Any,<:Simplex{2},<:Union{InverseDistanceCubed{2},InverseDistanceGradient{2}}},
    WeightedDomain{3,<:Any,<:SurfaceTriangle,<:Union{InverseDistanceCubed{3},InverseDistanceGradient{3}}}}

# The kernels are not invariant under affine maps (only under similarities), so the rule is
# built where it is asked for: such a domain is its own reference.
isreference(::PointKernelDomain) = true
reference(d::PointKernelDomain) = d

_exact(x::Integer) = Rational{BigInt}(x)
_exact(x::Rational) = Rational{BigInt}(x)
_exact(x::AbstractFloat) = Rational{BigInt}(x)      # exact: a float is a binary fraction
_exact(x::Real) = Rational{BigInt}(BigFloat(x))
_ulp(x::Real) = Rational{BigInt}(value_ulp(x))

# twice the signed area of (o, a, b) along the triangle's orientation `n` (in space), or in
# the plane; and twice the unsigned area of the parallelogram on a and b
_signed_area2(a::SVector{2}, b::SVector{2}, _) = a[1] * b[2] - a[2] * b[1]
_signed_area2(a::SVector{3}, b::SVector{3}, n) = dot(n, cross(a, b))
_area2(a::SVector{2}, b::SVector{2}) = abs(a[1] * b[2] - a[2] * b[1])
_area2(a::SVector{3}, b::SVector{3}) = norm(cross(a, b))

"""
    kernel_geometry(dom) -> (; x0, c, H2, patches, place, moved)

How a rule for `dom` is built, decided in exact arithmetic. `place` is where the singular
point lies: `:vertex`, `:edge` or `:interior` when it is on the triangle (within rounding, see
[`InverseDistance`](@ref)), and `:near` when it is off it. `x0` is the point the rule is for
(moved onto the triangle when it was within rounding of it, by `moved`; a point off the
triangle is taken as given). `c` is the apex the triangle is cut at: `x0` itself when it is on
the triangle, and otherwise the point of the triangle nearest to it, at the distance
`√H2 > 0`. `patches` are the sub-triangles `(c, p, q)` of positive area, as pairs `(p, q)` in
the triangle's orientation. Throws if the triangle is degenerate.
"""
function kernel_geometry(dom::PointKernelDomain)
    D = length(dom.weight.x0)
    V = [SVector{D,Rational{BigInt}}(map(_exact, p)) for p in vertices(dom.base)]
    given = SVector{D,Rational{BigInt}}(map(_exact, dom.weight.x0))
    tol2 = (64 * maximum(_ulp, vcat(collect(dom.weight.x0), (collect(p) for p in vertices(dom.base))...)))^2
    e1, e2 = V[2] - V[1], V[3] - V[1]
    n = D == 3 ? cross(e1, e2) : nothing
    total = _signed_area2(e1, e2, n)
    iszero(total) && throw(ArgumentError("the triangle $(dom.base) is degenerate"))
    edges = ((2, 3), (3, 1), (1, 2))
    λ(p) = [_signed_area2(V[j] - p, V[k] - p, n) / total for (j, k) in edges]
    # x projected onto the plane, and whether it was within rounding of it
    x, inplane = given, true
    if D == 3
        s = dot(n, given - V[1])
        nn = dot(n, n)
        x = given - (s / nn) * n
        inplane = s^2 <= tol2 * nn
    end
    # the nearest point of the boundary; a point within rounding of it is moved onto it, and
    # onto a vertex if as close to one
    c, δ2 = nothing, nothing
    for (j, k) in edges
        d = V[k] - V[j]
        p = V[j] + clamp(dot(x - V[j], d) / dot(d, d), 0, 1) * d
        q = dot(x - p, x - p)
        (δ2 === nothing || q < δ2) && ((c, δ2) = (p, q))
    end
    inside = all(>=(0), λ(x))
    if δ2 <= tol2
        x = c
        for v in V
            dot(x - v, x - v) <= tol2 && (x = v)
        end
    elseif !inside
        x = c                                       # the nearest point of the triangle
    end
    patches = [(V[j], V[k]) for ((j, k), l) in zip(edges, λ(x)) if l > 0]
    if inplane && (inside || δ2 <= tol2)            # on the triangle: singular
        place = length(patches) == 1 ? :vertex : length(patches) == 2 ? :edge : :interior
        return (; x0 = x, c = x, H2 = zero(Rational{BigInt}), patches, place,
                moved = sqrt(Float64(dot(x - given, x - given))))
    end
    return (; x0 = given, c = x, H2 = dot(given - x, given - x), patches, place = :near, moved = 0.0)
end

"The number of sub-triangles the triangle is cut into at the apex: 1, 2 or 3."
npatches(dom::PointKernelDomain) = length(kernel_geometry(dom).patches)

# On the sub-triangle (x₀, p, q) in polar coordinates about x₀ the kernel cancels the
# Jacobian, ∫∫ f (1/r) r dr dψ, and along the edge pq, at the foot f of the perpendicular from
# x₀ at distance h, the angle is ψ = atan(τ/h), dψ = h dτ/(h² + τ²) and the ray has length
# R = √(h² + τ²). So the patch contributes ∫ dτ h/√(h² + τ²) ∫₀¹ f(x₀ + ρ (e(τ) − x₀)) dρ;
# with f = 1 that is h [asinh(τ/h)] between the ends of the edge.
"The foot, unit direction, distance `h` and edge parameters `τ_p, τ_q` of a patch, in `S`."
function patch_frame(x0, p, q, ::Type{S}) where {S}
    D = length(x0)
    x0, p, q = SVector{D,S}(x0), SVector{D,S}(p), SVector{D,S}(q)
    u = (q - p) / norm(q - p)
    f = p + dot(x0 - p, u) * u
    return f, u, norm(x0 - f), dot(p - f, u), dot(q - f, u)
end

"""
`∫_T 1/|y − x₀| dy` in the ambient `BigFloat` precision: in closed form for a point on the
triangle, and in polar coordinates about the nearest point otherwise (see
`polar_kernel_integrals`).
"""
function measure(d::InverseDistanceDomain)
    g = kernel_geometry(d)
    g.place === :near && return only(polar_kernel_integrals(d, 1, BigFloat, (out, y) -> (out[1] = one(BigFloat))))
    return sum(g.patches) do (p, q)
        _, _, h, τp, τq = patch_frame(g.x0, p, q, BigFloat)
        h * (asinh(τq / h) - asinh(τp / h))
    end
end

"""
The finite part of `∫_T K(y) dy` for the kernel of a [`FinitePartDomain`](@ref), in the ambient
`BigFloat` precision (see `polar_finitepart_integrals`).
"""
measure(d::FinitePartDomain) = only(polar_finitepart_integrals(d, 1, BigFloat, (out, y) -> (out[1] = one(BigFloat))))
