# Applying rules (PLAN §3.1).
#
# There is deliberately no `integrate(f, domain; degree = 20)`: construction must be a
# visible subexpression, `integrate(f, rule(domain; degree = 20))`, so that it is hoisted
# out of loops. The package does not cache rules (§7), so such an overload would silently
# rebuild the rule on every call.
#
# Accumulation is generic in the integrand's return type: the accumulator starts from the
# first term, w₁ f(x₁), so its type is exactly `typeof(w[1] * f(x[1]))` (complex, vector,
# StaticArray, Unitful, …) and no `zero(T)` of the wrong type is ever formed. Nothing
# allocates when the nodes are isbits and `f` does not. The integrand is always typed
# `f::F` so that it is specialised on even though it is only passed through (Julia does not
# specialise on pass-through function arguments by default, which would force dynamic
# dispatch and box the node tuples).

"""
    integrate(f, rule; batch = false)
    integrate(f, rule, domain)
    integrate(f, rule, cells::AbstractVector{<:Domain})

`Σᵢ wᵢ f(xᵢ)`. In 1D `f` receives a scalar, otherwise an `SVector`.

- With a `domain`, the rule is mapped affinely onto it on the fly (no new rule is built).
- With a vector of `cells` (a mesh), the sum over all cells, rule reused.
- With `batch = true`, `f` is called once with all nodes — a `D × N` matrix (a vector in
  1D) — and must return the `N` values.

For hot loops use `static(rule)`, which keeps everything in registers.
"""
function integrate(f::F, r::QuadratureRule; batch::Bool = false) where {F}
    batch && return _integrate_batch(f, r)
    x = nodes(r)
    w = weights(r)
    @inbounds acc = w[1] * f(x[1])
    @inbounds for i in 2:length(w)
        acc += w[i] * f(x[i])
    end
    return acc
end

function _integrate_batch(f, r::AnyRule)
    x = nodes(r)
    X = eltype(x) <: Number ? collect(x) : reduce(hcat, x)
    vals = f(X)
    length(vals) == length(x) ||
        throw(DimensionMismatch("batched integrand returned $(length(vals)) values for $(length(x)) nodes"))
    w = weights(r)
    acc = w[1] * vals[1]
    for i in 2:length(w)
        acc += w[i] * vals[i]
    end
    return acc
end

function integrate(f::F, r::QuadratureRule, dom::Domain) where {F}
    A, b, J = _affine(domain(r), dom, _maptype(eltype(r), _coordtype(dom)))
    x = nodes(r)
    w = weights(r)
    @inbounds acc = (w[1] * J) * f(A * x[1] + b)
    @inbounds for i in 2:length(w)
        acc += (w[i] * J) * f(A * x[i] + b)
    end
    return acc
end

function integrate(f::F, r::AnyRule, cells::AbstractVector{<:Domain}; threaded::Bool = false) where {F}
    isempty(cells) && throw(ArgumentError("no cells to integrate over"))
    threaded && Threads.nthreads() > 1 && return _integrate_threaded(f, r, cells)
    acc = integrate(f, r, cells[1])
    for k in 2:length(cells)
        acc += integrate(f, r, cells[k])
    end
    return acc
end

# One accumulator per thread, summed in thread order: the result does not depend on how the
# cells were scheduled, only on how many threads there are.
#
# Every name assigned inside the threaded loop must be used nowhere else in this function.
# A name that is also assigned outside becomes one shared (boxed) variable, which every
# thread then writes — a race that silently returns wrong sums.
function _integrate_threaded(f::F, r::AnyRule, cells) where {F}
    nt = Threads.nthreads()
    partials = Vector{Any}(nothing, nt)
    chunks = [(t - 1) * length(cells) ÷ nt + 1:t * length(cells) ÷ nt for t in 1:nt]
    Threads.@threads :static for t in 1:nt
        rng = chunks[t]
        if !isempty(rng)
            part = integrate(f, r, cells[first(rng)])
            for k in (first(rng) + 1):last(rng)
                part += integrate(f, r, cells[k])
            end
            partials[t] = part
        end
    end
    done = [p for p in partials if p !== nothing]
    total = first(done)
    for j in 2:length(done)
        total += done[j]
    end
    return total
end

_coordtype(s::Simplex) = eltype(eltype(s.vertices))

"Arithmetic type of a mapped evaluation: exact inputs stay exact (Int → Rational{BigInt})."
@inline function _maptype(::Type{R}, ::Type{C}) where {R,C}
    S = promote_type(R, C)
    return S <: Integer ? Rational{BigInt} : S
end
_coordtype(d::Interval) = typeof(d.a)
_coordtype(d::Orthotope) = eltype(d.lo)

# ---------------------------------------------------------------------------------------
# The runtime form: inlined, and unrolled at compile time (the node count is a type
# parameter), so the node tuple is never indexed dynamically and never leaves registers.

@inline function integrate(f::F, r::StaticQuadratureRule; batch::Bool = false) where {F}
    batch && return _integrate_batch(f, r)
    return _unrolled_sum(f, r.nodes, r.weights)
end

@inline function integrate(f::F, r::StaticQuadratureRule, dom::Domain) where {F}
    A, b, J = _affine(domain(r), dom, _maptype(eltype(r), _coordtype(dom)))
    return J * _unrolled_sum(y -> f(A * y + b), r.nodes, r.weights)
end

@generated function _unrolled_sum(f::F, x::SVector{N}, w::SVector{N}) where {F,N}
    N > 256 && return :(_loop_sum(f, x, w))
    ex = :(w[1] * f(x[1]))
    for i in 2:N
        ex = :($ex + w[$i] * f(x[$i]))
    end
    return :(@inbounds $ex)
end

function _loop_sum(f::F, x, w) where {F}
    @inbounds acc = w[1] * f(x[1])
    @inbounds for i in 2:length(w)
        acc += w[i] * f(x[i])
    end
    return acc
end
