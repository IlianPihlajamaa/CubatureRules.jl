# Number types exercised in CI (PLAN §7): Float32, Float64, Double64, BigFloat at several
# precisions, Rational{BigInt} for the exact families.
using CubatureRules, Test
import CubatureRules: monomial_moment, ⊗   # public, not exported
import DoubleFloats: Double64   # `using` would make ⊗ ambiguous in Main
const CR = CubatureRules

@testset "Float32" begin
    r = rule(Simplex{2}(); degree = 8, T = Float32)
    @test eltype(r) === Float32
    @test passed(check(r))
    @test passed(check(rule(Interval(); degree = 9, T = Float32)))
end

@testset "Double64" begin
    r = rule(Simplex{2}(); degree = 12, T = Double64)
    @test eltype(r) === Double64
    v = check(r)
    @test v.exact && v.sharp === true
    @test integrate(x -> x[1]^5 * x[2]^6, r) ≈ Double64(factorial(5) * factorial(6) // factorial(big(13))) rtol = 1e-30
end

@testset "BigFloat at several precisions" begin
    for digits in (20, 50, 120)
        r = rule(Simplex{2}(); degree = 9, digits)
        @test precision(first(weights(r))) == CR.digits_to_bits(digits)
        @test passed(check(r))
        c = rule(ConicalProduct(), Simplex{2}(); degree = 9, digits)
        @test passed(check(c))
        g = rule(Interval(); degree = 9, digits)
        @test passed(check(g))
    end
end

@testset "Rational{BigInt} (exact families only)" begin
    r = rule(Simplex{3}(); degree = 3, T = Rational{BigInt})
    @test eltype(r) == Rational{BigInt}
    @test family(r) == "GrundmannMöller"
    @test integrate(x -> x[1]^2 * x[2], r) == monomial_moment(Simplex{3}(), (2, 1, 0))
    l = rule(Interval(); degree = 3, T = Rational{BigInt})
    @test family(l) == "NewtonCotes"          # the exact-rational 1D family
    @test eltype(l) == Rational{BigInt}
    @test sum(weights(l)) == 2
    @test integrate(x -> x^2, l) == 2 // 3
    @test_throws NoRuleError rule(WeightedDomain(Interval(), JacobiWeight(1, 0)); degree = 3, T = Rational{BigInt})
end
