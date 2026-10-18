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

# Over every kind of domain: the rule `rule` picks verifies; mapped onto an affine image of
# its domain it still does; and, chosen without `copyleft`, it carries no copyleft terms.
const DOMAINS = [Interval(), Simplex{2}(), Simplex{3}(), Orthotope{2}(), Orthotope{3}(), Disk(), Ball{3}(),
                 Sphere{3}(), Wedge(), Pyramid(), GaussianSpace(2), LaguerreRay(0), HermiteLine(),
                 WeightedDomain(Interval(), JacobiWeight(1, 2))]
domgen = Data.Integers(1, length(DOMAINS))      # an index into DOMAINS

@check max_examples = 40 function every_domain_verifies(i = domgen, d = Data.Integers(0, 9), big = Data.Booleans())
    dom = DOMAINS[i]
    r = big ? rule(dom; degree = d, digits = 30) : rule(dom; degree = d)
    # chosen without `copyleft`, it carries no copyleft terms (Lebedev.jl is loaded by now)
    degree(r) >= d && passed(check(r)) && !occursin("GPL", provenance(r).license)
end

# an affine image of the domain, with a positive Jacobian, from a few small integers
image(::Interval, a, b) = Interval(a, a + b)
image(::Simplex{2}, a, b) = Simplex((a, 0), (a + b, 1), (a, b + 1))
image(::Simplex{3}, a, b) = Simplex((a, 0, 0), (a + b, 0, 0), (a, b, 0), (a, 1, b))   # volume b³/6
image(::Orthotope{D}, a, b) where {D} = Orthotope(ntuple(i -> a + i, D), ntuple(i -> a + i + b, D))
image(::Ball{D}, a, b) where {D} = Ball(ntuple(i -> a - i, D), b)
image(::Sphere{D}, a, b) where {D} = Sphere(ntuple(i -> a + i, D), b)
image(dom, a, b) = nothing

@check max_examples = 40 function mapped_rules_verify(i = domgen, d = Data.Integers(0, 8), a = Data.Integers(-3, 3),
                                                       b = Data.Integers(1, 4))
    img = image(DOMAINS[i], a, b)
    img === nothing && return true
    r = rule(img; degree = d)
    passed(check(r)) && numbertype(r) === Float64
end

@check max_examples = 20 function subdivided_rules_verify(d = Data.Integers(0, 8), n = Data.Integers(1, 4),
                                                           tri = Data.Booleans())
    r = subdivide(rule(tri ? Simplex{2}() : Interval(); degree = d), n)
    passed(check(r))
end
