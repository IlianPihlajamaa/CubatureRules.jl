# Property-based tests over (family, degree, T) with Supposition.jl (PLAN §8): every rule
# the selector can build verifies against its claim, at family and precision boundaries.
using CubatureRules, Test, Supposition
import CubatureRules: ⊗   # public, not exported
import DoubleFloats: Double64   # `using` would make ⊗ ambiguous in Main
using Supposition: Data

const FAMILIES_2D = [XiaoGimbutas(), GrundmannMöller(), ConicalProduct()]
const TYPES = [Float32, Float64, Double64, BigFloat]

famgen = Data.SampledFrom(FAMILIES_2D)
typegen = Data.SampledFrom(TYPES)
deggen = Data.Integers(0, 20)

@check max_examples = 40 function rules_verify(f = famgen, d = deggen, T = typegen)
    r = T === BigFloat ? rule(f, Simplex{2}(); degree = d, digits = 40) : rule(f, Simplex{2}(); degree = d, T)
    v = check(r)
    degree(r) >= d && v.exact && v.sharp !== false && v.weights_sum_ok && v.symmetric !== false
end

@check max_examples = 20 function exact_rules_are_exact(D = Data.Integers(2, 4), s = Data.Integers(0, 3))
    r = rule(GrundmannMöller(), Simplex{D}(); degree = 2s + 1, T = Rational{BigInt})
    v = check(r)
    v.exact && v.max_residual == 0 && sum(weights(r)) == 1 // factorial(D)
end

@check max_examples = 20 function affine_maps_preserve_exactness(d = Data.Integers(1, 12),
                                                                  a = Data.Floats{Float64}(; infs = false, nans = false))
    x0 = clamp(a, -10.0, 10.0)
    t = Simplex((x0, 0.0), (x0 + 2.0, 0.5), (x0 + 0.3, 1.7))
    r = rule(t; degree = d)
    check(r).exact
end
