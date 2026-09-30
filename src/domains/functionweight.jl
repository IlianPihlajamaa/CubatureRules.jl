# A measure given by a function that can be evaluated at any precision (PLAN §6 Tier 2:
# Gautschi's multiple-component discretisation, `mcdis`).
#
# The measure is a sum of pieces and point masses. Each piece is a smooth function `g` times
# a classical factor that carries whatever makes the piece hard to integrate: an endpoint
# singularity (x − a)^β, exponential decay on a half line, or a Gaussian on the whole line.
# A piece is discretised by the Gauss rule of its classical factor, which integrates g times
# that factor to high accuracy with few points when g is smooth, and the discretised pieces
# and the point masses together form one discrete measure. Its recurrence coefficients come
# from the discrete Stieltjes procedure (families/onedim/stieltjes.jl), and the number of
# points is doubled until two resolutions agree to the precision asked for.

"""
    WeightPiece

One piece of a [`FunctionWeight`](@ref): `g(x)` times a classical factor, `kind` being
`:jacobi` (`(b − x)^α (x − a)^β` on `[a, b]`), `:laguerre` (`(x − a)^β e^{−λ(x − a)}` on
`[a, ∞)`) or `:hermite` (`e^{−x²}` on the real line).
"""
struct WeightPiece
    g::Any                     # unspecialised, for the reason given at MonicRecurrence
    kind::Symbol
    a::Real
    b::Real
    α::Real
    β::Real
    rate::Real
end

"""
    FunctionWeight(g, a, b; α = 0, β = 0, rate = 1, label)

A weight given by a function. One piece, on an interval chosen by `a` and `b`:

| `a`, `b` | weight |
|---|---|
| both finite | `g(x) (b − x)^α (x − a)^β` on `[a, b]` |
| `b = Inf` | `g(x) (x − a)^β e^{−rate (x − a)}` on `[a, ∞)` |
| `a = -Inf`, `b = Inf` | `g(x) e^{−x²}` on the real line |

Pieces and point masses combine with `+`, and [`PointMass`](@ref) adds a point mass:

```julia
w = FunctionWeight(x -> exp(-x), 0, 1; β = -1/2) + FunctionWeight(x -> one(x), 1, 2) + PointMass(2, 1/10)
rule(WeightedDomain(Interval(0, 2), w); degree = 19, digits = 40)
```

Each piece is discretised with the Gauss rule of its own factor, so put any endpoint
singularity or decay into `α`, `β` or the half-line factor and keep `g` smooth: then few
points suffice. A kink or a singularity left inside `g` makes the discretisation converge
slowly, and the construction reports that rather than returning an inaccurate rule. Split
the piece at the kink instead.

`g` is evaluated at `BigFloat` arguments at the working precision and must return a
`BigFloat` there: a function that returns `Float64` describes the weight only to `Float64`,
and is refused.
"""
struct FunctionWeight
    pieces::Vector{WeightPiece}
    atoms::Vector{Tuple{Real,Real}}   # (location, mass)
    label::String
end

function FunctionWeight(g, a::Real, b::Real; α::Real = 0, β::Real = 0, rate::Real = 1, label::AbstractString = "")
    a < b || throw(ArgumentError("a piece needs a < b, got [$a, $b]"))
    # a named function says what the weight is; an anonymous one can only be called g
    name = g isa Function && !startswith(string(nameof(g)), "#") ? string(nameof(g)) : "g"
    if isfinite(a) && isfinite(b)
        (α > -1 && β > -1) || throw(ArgumentError("the exponents must exceed −1 for the piece to be integrable"))
        piece = WeightPiece(g, :jacobi, a, b, α, β, rate)
        default = "$name(x)" * _factor(" (%s − x)^%s", b, α) * _factor(" (x − %s)^%s", a, β) * " on [$a, $b]"
    elseif isfinite(a) && b == Inf
        α == 0 || throw(ArgumentError("on a half line the exponent belongs to the finite end: use β"))
        (β > -1 && rate > 0) || throw(ArgumentError("need β > −1 and rate > 0"))
        piece = WeightPiece(g, :laguerre, a, b, α, β, rate)
        default = "$name(x)" * _factor(" (x − %s)^%s", a, β) * " e^{−$(rate == 1 ? "" : "$rate ")(x − $a)} on [$a, ∞)"
    elseif a == -Inf && b == Inf
        (α == 0 && β == 0) || throw(ArgumentError("a piece on the real line takes no exponents"))
        piece = WeightPiece(g, :hermite, a, b, α, β, rate)
        default = "$name(x) e^{−x²} on ℝ"
    else
        throw(ArgumentError("a piece on (−∞, $b] is not supported; substitute x → −x"))
    end
    return FunctionWeight([piece], Tuple{Real,Real}[], isempty(label) ? default : String(label))
end

_factor(fmt, e, p) = iszero(p) ? "" : replace(replace(fmt, "%s" => string(e); count = 1), "%s" => string(p); count = 1)

"""
    PointMass(x, m)

A point mass `m` at `x`, to add to a [`FunctionWeight`](@ref) with `+`.
"""
function PointMass(x::Real, m::Real)
    (isfinite(x) && m > 0) || throw(ArgumentError("a point mass needs a finite location and a positive mass"))
    return FunctionWeight(WeightPiece[], [(x, m)], "$m δ(x − $x)")
end

Base.:+(u::FunctionWeight, v::FunctionWeight) =
    FunctionWeight(vcat(u.pieces, v.pieces), vcat(u.atoms, v.atoms), u.label * " + " * v.label)
Base.show(io::IO, w::FunctionWeight) = print(io, w.label)

"The closed hull of the weight's support, as `(lo, hi)`; unbounded ends are `±Inf`."
function support(w::FunctionWeight)
    lo = minimum(vcat([p.a for p in w.pieces], [x for (x, _) in w.atoms]))
    hi = maximum(vcat([p.b for p in w.pieces], [x for (x, _) in w.atoms]))
    return lo, hi
end

"""
    FunctionDomain

A one-dimensional domain carrying a [`FunctionWeight`](@ref).
"""
const FunctionDomain = WeightedDomain{1,<:Any,<:Domain{1},<:FunctionWeight}

# Like a moment weight, a function weight is stated in absolute coordinates: the domain is its
# own reference and is never mapped.
isreference(::FunctionDomain) = true
reference(d::FunctionDomain) = d

"""
    check_support(dom::FunctionDomain)

Refuse a pairing in which the weight lives outside the domain it was paired with.
"""
function check_support(dom::FunctionDomain)
    lo, hi = support(dom.weight)
    blo, bhi = endpoints(dom.base)
    (blo <= lo && hi <= bhi) || throw(ArgumentError(
        "$(dom.weight) lives on [$lo, $hi], outside $(dom.base); pair it with a domain that contains it"))
    return nothing
end

"""
    discretize(w::FunctionWeight, M, bits) -> (x, λ)

The measure as a discrete one, in BigFloat at `bits`: an `M`-point Gauss rule of each piece's
classical factor with its weights multiplied by `g`, followed by the point masses.
"""
function discretize(w::FunctionWeight, M::Integer, bits::Integer)
    xs, λs = BigFloat[], BigFloat[]
    with_bits(bits) do
        for p in w.pieces
            if p.kind === :jacobi
                t, v, _ = gauss_jacobi_work(M, p.α, p.β, bits)      # (1 − t)^α (1 + t)^β on [−1, 1]
                a, b = BigFloat(p.a), BigFloat(p.b)
                h = (b - a) / 2
                scale = h^(BigFloat(p.α) + BigFloat(p.β) + 1)
                x = [a + h * (1 + ti) for ti in t]
                λ = [scale * vi for vi in v]
            elseif p.kind === :laguerre
                t, v, _ = gauss_from_recurrence(M, laguerre_recurrence(p.β), bits)
                r = BigFloat(p.rate)
                x = [BigFloat(p.a) + ti / r for ti in t]
                λ = [vi / r^(BigFloat(p.β) + 1) for vi in v]
            else
                x, λ, _ = gauss_from_recurrence(M, hermite_recurrence(), bits)
            end
            for i in eachindex(x)
                gx = p.g(x[i])
                gx isa BigFloat || throw(ArgumentError(
                    "the weight function returned a $(typeof(gx)) for a BigFloat argument, so it describes " *
                    "the weight only to that precision; write it so that it computes in the type of its argument"))
                push!(xs, x[i])
                push!(λs, λ[i] * gx)
            end
        end
        for (x, m) in w.atoms
            push!(xs, BigFloat(x))
            push!(λs, BigFloat(m))
        end
    end
    return xs, λs
end

"""
    discretized_moments(w, aux, K, bits) -> Vector{BigFloat}

The modified moments `∫ πₖ dλ`, `k = 0:K`, of `w` against the monic family `aux`, by
discretisations with increasing numbers of points until two agree to `bits` bits. For
verification: independent of the Stieltjes procedure that builds the rules.
"""
function discretized_moments(w::FunctionWeight, aux::MonicRecurrence, K::Integer, bits::Integer;
                             maxpoints::Integer = 4096)
    function moments(M)
        x, λ = discretize(w, M, bits + 32)
        with_bits(bits + 32) do
            m, scale = zeros(BigFloat, K + 1), zeros(BigFloat, K + 1)
            for (xi, λi) in zip(x, λ)
                p0, p1 = zero(BigFloat), one(BigFloat)
                for k in 0:K
                    m[k + 1] += λi * p1
                    scale[k + 1] += abs(λi * p1)
                    p0, p1 = p1, (xi - aux.a(k, BigFloat)) * p1 - aux.b(k, BigFloat) * p0
                end
            end
            m, scale
        end
    end
    M = max(2K + 2, 16)
    prev, _ = moments(M)
    while 2M <= maxpoints
        M *= 2
        cur, scale = moments(M)
        all(k -> abs(cur[k] - prev[k]) <= ldexp(scale[k], -(bits + 4)), eachindex(cur)) && return cur
        prev = cur
    end
    throw(RefinementError("StieltjesDiscretization",
        "the moments of $(w) did not converge with $maxpoints points per piece; is g smooth on each piece?"))
end

"""
    verification_weight(w::FunctionWeight) -> MomentWeight

`w` described by its modified moments against monic polynomials suited to its support
(Legendre on a finite hull, Laguerre on a half line, Hermite on the line), each moment
computed by discretisation. [`verify`](@ref) tests rules on a [`FunctionDomain`](@ref)
against these.
"""
function verification_weight(w::FunctionWeight)
    lo, hi = support(w)
    aux = if isfinite(lo) && isfinite(hi)
        shift(monic(jacobi_recurrence(0, 0)), lo, hi)
    elseif isfinite(lo)
        MonicRecurrence((k, T) -> T(2k + 1) + T(lo), (k, T) -> T(k)^2)
    else
        MonicRecurrence((k, T) -> zero(T), (k, T) -> T(k) / 2)
    end
    cache = Dict{Tuple{Int,Int},Vector{BigFloat}}()
    lock = ReentrantLock()
    function moment(k, ::Type{T}) where {T}
        prec = precision(BigFloat)
        m = Base.lock(lock) do
            key = (prec, 16 * cld(k + 1, 16))
            get!(() -> discretized_moments(w, aux, key[2], prec), cache, key)
        end
        return T(m[k + 1])
    end
    return MomentWeight(aux, moment; label = w.label)
end

measure(d::FunctionDomain) = verification_weight(d.weight).moments(0, BigFloat)
