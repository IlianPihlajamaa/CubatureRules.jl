using CubatureRules, Test
const CR = CubatureRules

# Weights given as functions: Gautschi's multiple-component discretisation with the discrete
# Stieltjes procedure. Classical weights written as function weights must reproduce the
# classical Gauss rules exactly; a new weight must agree with the independent modified
# Chebyshev algorithm run on its exact moments.

samerule(r, s) = nodes(r) == nodes(s) && weights(r) == weights(s)
unit(x) = one(x)

@testset "function weights reproduce the classical rules" begin
    d = 30
    @test samerule(rule(WeightedDomain(Interval(-1, 1), FunctionWeight(unit, -1, 1)); degree = 11, digits = d),
                   rule(GaussLegendre(), Interval(); degree = 11, digits = d))
    @test samerule(rule(WeightedDomain(Interval(-1, 1), FunctionWeight(unit, -1, 1; α = 1 // 2, β = -1 // 3)); degree = 11, digits = d),
                   rule(GaussJacobi(1 // 2, -1 // 3), WeightedDomain(Interval(), JacobiWeight(1 // 2, -1 // 3)); degree = 11, digits = d))
    # two pieces that together make the Legendre weight
    split = FunctionWeight(unit, -1, 0) + FunctionWeight(unit, 0, 1)
    @test samerule(rule(WeightedDomain(Interval(-1, 1), split); degree = 11, digits = d),
                   rule(GaussLegendre(), Interval(); degree = 11, digits = d))
    # half line and real line
    lag = rule(WeightedDomain(HalfLine(), FunctionWeight(unit, 0, Inf; β = 1 // 2)); degree = 11, digits = d)
    @test samerule(lag, rule(GaussLaguerre(1 // 2), LaguerreRay(1 // 2); degree = 11, digits = d))
    @test samerule(rule(WeightedDomain(RealLine(), FunctionWeight(unit, -Inf, Inf)); degree = 11, digits = d),
                   rule(GaussHermite(), HermiteLine(); degree = 11, digits = d))
    # e^{−2x}: the Laguerre rule scaled by 1/2
    lag2 = rule(WeightedDomain(HalfLine(), FunctionWeight(unit, 0, Inf; rate = 2)); degree = 11, digits = d)
    ref = rule(GaussLaguerre(), LaguerreRay(); degree = 11, digits = d)
    @test maximum(abs, 2 .* nodes(lag2) .- nodes(ref)) < 1e-28
    @test maximum(abs, 2 .* weights(lag2) .- weights(ref)) < 1e-28
end

@testset "a new weight agrees with modified Chebyshev on its exact moments" begin
    # e^{−x} on [0, 2] plus a point mass 1/2 at 1: ordinary moments γ(k + 1, 2) + 1/2
    w = FunctionWeight(x -> exp(-x), 0, 2) + PointMass(1, 1 // 2)
    r = rule(WeightedDomain(Interval(0, 2), w); degree = 15, digits = 40)
    @test family(r) == "StieltjesDiscretization" && npoints(r) == 8
    lowergamma(k, T) = factorial(big(k)) * (1 - exp(-T(2)) * sum(T(2)^j / factorial(big(j)) for j in 0:k))
    mw = OrdinaryMoments((k, T) -> lowergamma(k, T) + T(1 // 2); support = (0, 2))
    rm = rule(WeightedDomain(Interval(0, 2), mw); degree = 15, digits = 40)
    @test maximum(abs, nodes(r) .- nodes(rm)) < 1e-38
    @test maximum(abs, weights(r) .- weights(rm)) < 1e-38
    @test passed(verify(r))
    @test occursin("Gautschi", sprint(show, MIME"text/plain"(), r))
    # and in Float64
    r64 = rule(WeightedDomain(Interval(0, 2), w); degree = 11)
    @test eltype(r64) == Float64 && passed(verify(r64))
    @test measure(WeightedDomain(Interval(0, 2), w)) ≈ 1 - exp(-2) + 1 / 2
end

@testset "singularities belong in the exponents" begin
    # √x as an exponent: fast and verified
    r = rule(WeightedDomain(Interval(0, 1), FunctionWeight(unit, 0, 1; β = 1 // 2)); degree = 9, digits = 30)
    @test passed(verify(r))
    # √x inside g converges only algebraically, and the construction says so quickly
    t = @elapsed err = try
        rule(WeightedDomain(Interval(0, 1), FunctionWeight(sqrt, 0, 1)); degree = 9, digits = 30)
        nothing
    catch e
        e
    end
    @test err isa RefinementError && occursin("exponents", sprint(showerror, err))
    @test t < 30
end

@testset "what a function weight refuses" begin
    @test_throws ArgumentError rule(WeightedDomain(Interval(0, 1), FunctionWeight(x -> sin(Float64(x)), 0, 1)); degree = 5, digits = 30)
    @test_throws ArgumentError rule(WeightedDomain(Interval(0, 1), FunctionWeight(exp, 0, 2)); degree = 5)
    @test_throws ArgumentError FunctionWeight(unit, 1, 0)
    @test_throws ArgumentError FunctionWeight(unit, 0, 1; β = -1)
    @test_throws ArgumentError FunctionWeight(unit, -Inf, 0)
    @test_throws ArgumentError PointMass(0, -1)
    # point masses alone: rules with fewer points than masses exist, others are refused
    atoms = PointMass(0.25, 1) + PointMass(0.5, 2) + PointMass(0.75, 1)
    r = rule(StieltjesDiscretization(), WeightedDomain(Interval(0, 1), atoms); npoints = 2, digits = 30)
    @test integrate(x -> x^3, r) ≈ 0.25^3 + 2 * 0.5^3 + 0.75^3
    @test_throws NoRuleError rule(StieltjesDiscretization(), WeightedDomain(Interval(0, 1), atoms); npoints = 3)
    @test string(FunctionWeight(exp, 0, 1; β = 1 // 2)) == "exp(x) (x − 0)^1//2 on [0, 1]"
end
