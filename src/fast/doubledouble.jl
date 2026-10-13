# Double-double arithmetic for the few quantities of the fast kernel rules whose rounding in
# Float64 would be amplified (notes/fast-kernel-rules.md): sums and products by error-free
# transformations, `asinh` and `log` through MultiFloats' Float64x2. A value is a pair
# (hi, lo) with |lo| ≤ ulp(hi)/2.

@inline function two_sum(a::Float64, b::Float64)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end
@inline two_prod(a::Float64, b::Float64) = (p = a * b; (p, fma(a, b, -p)))

@inline function dd_add((ah, al)::NTuple{2,Float64}, (bh, bl)::NTuple{2,Float64})
    s, e = two_sum(ah, bh)
    return two_sum(s, e + al + bl)
end
@inline dd_neg((h, l)::NTuple{2,Float64}) = (-h, -l)
@inline function dd_mul((ah, al)::NTuple{2,Float64}, (bh, bl)::NTuple{2,Float64})
    p, e = two_prod(ah, bh)
    return two_sum(p, e + ah * bl + al * bh)
end
"`(ah + al)/b` for a Float64 `b`."
@inline function dd_div((ah, al)::NTuple{2,Float64}, b::Float64)
    q = ah / b
    r = fma(-q, b, ah) + al
    return two_sum(q, r / b)
end

_x2((h, l)::NTuple{2,Float64}) = Float64x2(h) + Float64x2(l)
function _dd(y::Float64x2)
    h = Float64(y)
    return (h, Float64(y - Float64x2(h)))
end
function dd_asinh(x::NTuple{2,Float64})
    y = _x2(x)
    a = abs(y)
    r = log(a + sqrt(a * a + 1))
    return _dd(y < 0 ? -r : r)
end
dd_log(x::NTuple{2,Float64}) = _dd(log(_x2(x)))
