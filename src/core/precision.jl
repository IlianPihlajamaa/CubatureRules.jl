# Precision requests and the build context (PLAN §2.7).
#
# Every generator computes in BigFloat at an explicitly passed precision (target + guard)
# and rounds exactly once at the end. MPFR is correctly rounded, which is what makes
# `rule(...)` bitwise reproducible across platforms.

"""
    CancellationToken()

Passed as `cancel = token` to [`rule`](@ref). The Newton driver checks it once per
iteration; [`cancel!`](@ref) makes the construction throw a [`CancelledError`](@ref) at
the next check.
"""
struct CancellationToken
    flag::Threads.Atomic{Bool}
end
CancellationToken() = CancellationToken(Threads.Atomic{Bool}(false))
cancel!(t::CancellationToken) = (t.flag[] = true; t)
iscancelled(t::CancellationToken) = t.flag[]
iscancelled(::Nothing) = false

"Thrown when a construction is cancelled through its [`CancellationToken`](@ref)."
struct CancelledError <: Exception end
Base.showerror(io::IO, ::CancelledError) = print(io, "rule construction was cancelled")

checkcancel(t) = iscancelled(t) && throw(CancelledError())

"""
    BuildContext

Everything a family's `build` method needs besides the family, domain and degree: the
output number type `T`, the target precision in bits, and the cancellation token.
"""
struct BuildContext{T}
    bits::Int
    cancel::Union{Nothing,CancellationToken}
    verbose::Int
end
BuildContext{T}(bits::Int; cancel = nothing, verbose = 0) where {T} =
    BuildContext{T}(bits, cancel, verbosity(verbose))

"""
    verbosity(v) -> Int

Normalise a `verbose` argument. `false` and `0` are silent; `true` and `1` report what is
about to be attempted, each iteration of an iterative solve, and anything that goes wrong;
`2` adds the inner detail — line-search backtracking, precision escalation, guard re-runs.

Progress goes through `@info`, so it obeys the ambient logger and can be captured, filtered
or redirected like any other Julia logging.
"""
verbosity(v::Bool) = v ? 1 : 0
verbosity(v::Integer) = Int(v)

outtype(::BuildContext{T}) where {T} = T
isexact(::BuildContext{T}) where {T} = T <: Rational || T <: Integer
target_digits(ctx::BuildContext) = isexact(ctx) ? typemax(Int) : floor(Int, ctx.bits * log10(2))

bits_of(::Type{Float16}) = 11
bits_of(::Type{Float32}) = 24
bits_of(::Type{Float64}) = 53
bits_of(::Type{BigFloat}) = precision(BigFloat)
bits_of(::Type{T}) where {T<:AbstractFloat} = try
    precision(T)
catch
    ceil(Int, -log2(eps(T))) + 1
end

digits_to_bits(d::Integer) = ceil(Int, d * log2(10)) + 1

"""
    resolve_precision(T, digits) -> (T, bits)

`digits` alone selects `BigFloat` at that many decimal digits; `T = BigFloat` alone uses
the ambient BigFloat precision *at the time of the call*; any other `T` uses its own
precision. Nothing downstream reads the ambient precision again.
"""
function resolve_precision(T, digits)
    if digits !== nothing
        T === nothing || T === BigFloat ||
            throw(ArgumentError("`digits` selects BigFloat; do not combine it with T = $T"))
        digits >= 1 || throw(ArgumentError("digits must be positive"))
        return BigFloat, digits_to_bits(digits)
    end
    T === nothing && return Float64, 53
    T <: Rational && return T, 0
    T <: AbstractFloat || throw(ArgumentError("unsupported number type $T"))
    return T, bits_of(T)
end

"""
    with_bits(f, bits)

Run `f()` with BigFloat precision set to `bits`. The only place the package sets
precision; generators receive `bits` explicitly and never read the ambient value.
"""
with_bits(f, bits::Integer) = setprecision(f, BigFloat, bits)

"Round a working-precision value once, to the output type of the context."
finalize_number(ctx::BuildContext{BigFloat}, x) = BigFloat(x; precision = ctx.bits)
finalize_number(ctx::BuildContext{T}, x) where {T} = T(x)
