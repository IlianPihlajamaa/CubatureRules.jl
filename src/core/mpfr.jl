# In-place MPFR arithmetic for the hot loops of the BigFloat bases and moment systems.
#
# Every BigFloat operation allocates its result and looks up the default precision, which
# since Julia 1.11 lives in a ScopedValue. In the Dubiner recurrences that overhead was about
# eight times the arithmetic, and it made BigFloat refinement at high degree take minutes.
# These wrappers write into an existing BigFloat instead.
#
# Two rules follow from that. The result is rounded to the precision of the *destination*,
# so every destination must be created at the working precision. And a destination must
# not be shared: `zeros(BigFloat, n)` fills a vector with one shared object, so buffers are
# made with `bigfloats`, which creates a distinct number for every entry.

const MPFR_RN = Base.MPFR.MPFRRoundNearest
const MPFR_LIB = Base.MPFR.libmpfr

for (f, c) in ((:mp_mul!, :mpfr_mul), (:mp_add!, :mpfr_add), (:mp_sub!, :mpfr_sub), (:mp_div!, :mpfr_div))
    @eval @inline function $f(z::BigFloat, x::BigFloat, y::BigFloat)
        ccall(($(QuoteNode(c)), MPFR_LIB), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
              z, x, y, MPFR_RN)
        return z
    end
end

"`z = x y + a`, rounded once. `z` may be `a`."
@inline function mp_fma!(z::BigFloat, x::BigFloat, y::BigFloat, a::BigFloat)
    ccall((:mpfr_fma, MPFR_LIB), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, y, a, MPFR_RN)
    return z
end

"`z = x y − a`, rounded once."
@inline function mp_fms!(z::BigFloat, x::BigFloat, y::BigFloat, a::BigFloat)
    ccall((:mpfr_fms, MPFR_LIB), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, y, a, MPFR_RN)
    return z
end

@inline function mp_set!(z::BigFloat, x::BigFloat)
    ccall((:mpfr_set, MPFR_LIB), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, MPFR_RN)
    return z
end

@inline function mp_set_si!(z::BigFloat, k::Integer)
    ccall((:mpfr_set_si, MPFR_LIB), Int32, (Ref{BigFloat}, Clong, Base.MPFR.MPFRRoundingMode), z, k, MPFR_RN)
    return z
end

"`z = k − x` for an integer `k`."
@inline function mp_si_sub!(z::BigFloat, k::Integer, x::BigFloat)
    ccall((:mpfr_si_sub, MPFR_LIB), Int32, (Ref{BigFloat}, Clong, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
          z, k, x, MPFR_RN)
    return z
end

"`z = x − k` for an integer `k`."
@inline function mp_sub_si!(z::BigFloat, x::BigFloat, k::Integer)
    ccall((:mpfr_sub_si, MPFR_LIB), Int32, (Ref{BigFloat}, Ref{BigFloat}, Clong, Base.MPFR.MPFRRoundingMode),
          z, x, k, MPFR_RN)
    return z
end

@inline function mp_sqrt!(z::BigFloat, x::BigFloat)
    ccall((:mpfr_sqrt, MPFR_LIB), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, MPFR_RN)
    return z
end

"`z = 2x`, exactly."
@inline function mp_twice!(z::BigFloat, x::BigFloat)
    ccall((:mpfr_mul_2ui, MPFR_LIB), Int32, (Ref{BigFloat}, Ref{BigFloat}, Culong, Base.MPFR.MPFRRoundingMode),
          z, x, 1, MPFR_RN)
    return z
end

"Distinct zeros of type `S`, one object per entry (for `BigFloat`, at the current precision)."
bigfloats(::Type{S}, dims::Integer...) where {S} = S === BigFloat ? [BigFloat(0) for _ in CartesianIndices(dims)] :
                                                   zeros(S, dims...)

"""
    unshare!(v, prec)

Make every entry of the `BigFloat` array `v` a distinct number of precision `prec`, so the
in-place kernels can write into it. Cheap once it holds: it then only compares.
"""
function unshare!(v::AbstractArray{BigFloat}, prec::Integer)
    prev = nothing                     # the previous entry as it was, not as replaced
    @inbounds for i in eachindex(v)
        cur = isassigned(v, i) ? v[i] : nothing
        if cur === nothing || cur === prev || precision(cur) != prec
            v[i] = BigFloat(0; precision = prec)
        end
        prev = cur
    end
    return v
end
