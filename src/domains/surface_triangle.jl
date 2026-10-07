# A flat triangle in three dimensions: a boundary element of a triangulated surface.

"""
    SurfaceTriangle(v₀, v₁, v₂)

A flat triangle in three-dimensional space, such as a boundary element of a triangulated
surface. A rule on it is a rule on the reference triangle `Simplex{2}()` carried over by the
affine map `ξ ↦ v₀ + ξ₁ (v₁ − v₀) + ξ₂ (v₂ − v₀)`, with the weights scaled by the area
element `|(v₁ − v₀) × (v₂ − v₀)|`; a degree claim means exactness for polynomials in the
three coordinates, on the triangle.

```julia
T = SurfaceTriangle((0, 0, 0), (1, 0, 1), (0, 1, 1))
r = rule(T; degree = 6)
integrate(y -> y[3]^2, r)
```

`integrate(f, r, T)` maps a rule `r` on `Simplex{2}()` onto `T` on the fly, without building a
new rule: the way to integrate over many elements with one rule. With an
[`InverseDistance`](@ref) kernel the triangle carries the weakly singular integrals of
boundary-element methods.
"""
struct SurfaceTriangle{T} <: Domain{3,T}
    vertices::SVector{3,SVector{3,T}}
    function SurfaceTriangle{T}(v::SVector{3,SVector{3,T}}) where {T}
        iszero(cross(v[2] - v[1], v[3] - v[1])) &&
            throw(ArgumentError("the triangle $(Tuple(Tuple.(v))) is degenerate: its vertices are collinear"))
        return new{T}(v)
    end
end
function SurfaceTriangle(v0, v1, v2)
    all(v -> length(v) == 3, (v0, v1, v2)) ||
        throw(ArgumentError("a triangle in three dimensions needs three coordinates per vertex"))
    T = promote_type(map(typeof, (v0..., v1..., v2...))...)
    return SurfaceTriangle{T}(SVector{3,SVector{3,T}}(SVector{3,T}(Tuple(v0)), SVector{3,T}(Tuple(v1)),
                                                      SVector{3,T}(Tuple(v2))))
end

vertices(t::SurfaceTriangle) = t.vertices
isreference(::SurfaceTriangle) = false
reference(::SurfaceTriangle) = Simplex{2}()
convert_domain(::Type{S}, t::SurfaceTriangle) where {S} = SurfaceTriangle{S}(map(v -> SVector{3,S}(map(S, v)), t.vertices))
Base.show(io::IO, t::SurfaceTriangle) = print(io, "SurfaceTriangle(", join((Tuple(v) for v in t.vertices), ", "), ")")

"""
    surface_frame(t, S) -> (E, v₀, J)

The edge matrix `E = [v₁ − v₀  v₂ − v₀]` (3 × 2), the vertex `v₀` and the area element
`J = |(v₁ − v₀) × (v₂ − v₀)|`, twice the area, in type `S`. `J` is a square root: in exact
arithmetic it exists only when it is rational, and otherwise this throws.
"""
function surface_frame(t::SurfaceTriangle, ::Type{S}) where {S}
    v0, v1, v2 = (SVector{3,S}(map(S, p)) for p in t.vertices)
    e1, e2 = v1 - v0, v2 - v0
    n = cross(e1, e2)
    return hcat(e1, e2), v0, _sqrt_of(dot(n, n), t)          # allocation-free for isbits S
end
_sqrt_of(x::AbstractFloat, _) = sqrt(x)
function _sqrt_of(x::Rational, t)
    a, b = isqrt(numerator(x)), isqrt(denominator(x))
    (a^2 == numerator(x) && b^2 == denominator(x)) ||
        throw(ArgumentError("the area of $t is irrational, so a rule on it cannot have exact rational weights; " *
                            "ask for floating-point output (T = Float64, or digits)"))
    return a // b
end

function measure(t::SurfaceTriangle{T}) where {T}
    S = T <: AbstractFloat ? T : BigFloat
    return surface_frame(t, S)[3] / 2
end

"Barycentric coordinates of the orthogonal projection of `x` onto the plane of `t`."
function barycentric(t::SurfaceTriangle, x::AbstractVector)
    S = float(promote_type(eltype(eltype(t.vertices)), eltype(x)))
    E, v0, _ = surface_frame(t, S)
    ξ = (E' * E) \ (E' * (SVector{3,S}(map(S, x)) - v0))
    return SVector{3,S}(1 - ξ[1] - ξ[2], ξ[1], ξ[2])
end

"""
    value_ulp(x)

The unit in the last place of `x` at the precision `x` carries: `eps(x)` for a machine float,
and for a `BigFloat` its own precision, not the ambient one that `eps` uses. Zero for exact
numbers (integers and rationals), which carry no rounding.
"""
value_ulp(x::AbstractFloat) = eps(x)
value_ulp(x::BigFloat) = iszero(x) || !isfinite(x) ? zero(BigFloat) : ldexp(BigFloat(1), exponent(x) - precision(x) + 1)
value_ulp(::Real) = 0

# A point computed in floating point is in the plane only to rounding: allow 64 units in the
# last place of its largest coordinate, at the precision the coordinates carry (a BigFloat
# node has its own precision, whatever the ambient one).
function _near_plane(t::SurfaceTriangle, x)
    S = float(promote_type(eltype(eltype(t.vertices)), eltype(x)))
    v = [SVector{3,S}(map(S, p)) for p in t.vertices]
    n = cross(v[2] - v[1], v[3] - v[1])
    ulp = maximum(c -> S(value_ulp(c)), vcat(collect(x), (collect(p) for p in t.vertices)...))
    return abs(dot(n, SVector{3,S}(map(S, x)) - v[1])) <= 64 * ulp * norm(n)
end
indomain(x, t::SurfaceTriangle; tol = 0) = _near_plane(t, x) && all(>=(-tol), barycentric(t, x))
isinterior(x, t::SurfaceTriangle; tol = 0) = _near_plane(t, x) && all(>(tol), barycentric(t, x))
