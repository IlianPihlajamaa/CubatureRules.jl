using CubatureRules, Test
const CR = CubatureRules

# Filon-type rules for e^{iωx}, cos ωx and sin ωx. The integrals of e^{βx} are known in
# closed form, (e^{zb} − e^{za})/z with z = β + iω, and their real and imaginary parts.

function osc_exact(kind, β, ω, a, b)
    setprecision(BigFloat, 400) do
        z = big(β) + im * big(ω)
        v = (exp(z * big(b)) - exp(z * big(a))) / z
        kind === :exp ? v : kind === :cos ? real(v) : imag(v)
    end
end
osc_weight(kind, ω) = kind === :exp ? Oscillatory(ω) : kind === :cos ? Oscillatory(ω, cos) : Oscillatory(ω, sin)
integral(r, β) = setprecision(() -> integrate(x -> exp(β * x), r), BigFloat, 400)

@testset "spherical Bessel functions both ways" begin
    setprecision(BigFloat, 256) do
        for s in ("1e-30", "1e-8", "0.3", "7.5", "-40", "1000", "123456.7"), N in (5, 60)
            κ = parse(BigFloat, s)
            up, down = CR.spherical_bessel_up(N, κ), CR.spherical_bessel_down(N, κ)
            @test all(abs(u - d) <= 1e-70 * max(abs(d), floatmin(Float64)) for (u, d) in zip(up, down))
        end
        # closed forms
        κ = big"2.5"
        j = CR.spherical_bessel_down(2, κ)
        @test abs(j[3] - ((3 / κ^2 - 1) * sin(κ) / κ - 3cos(κ) / κ^2)) < 1e-70
    end
end

@testset "Filon rules: closed forms, verification, output types" begin
    β = 7 // 10
    for kind in (:exp, :cos, :sin), (ω, a, b) in ((3, -1, 1), (50, 0, 2), (-200, 1, 3), (1 // 10^6, 0, 1))
        dom = WeightedDomain(Interval(a, b), osc_weight(kind, ω))
        r = rule(dom; degree = 40, digits = 30)
        @test family(r) == "Filon" && npoints(r) == 41 && degree(r) >= 40
        @test eltype(r) == BigFloat && (eltype(weights(r)) <: Complex) == (kind === :exp)
        @test first(nodes(r)) == a && last(nodes(r)) == b
        v = verify(r)
        @test passed(v) && v.exact && !v.interior
        @test abs(integral(r, β) - osc_exact(kind, β, ω, a, b)) < 1e-30
    end
    for T in (Float64, Float32), k in (Oscillatory(30), Oscillatory(30, cos), Oscillatory(30, sin))
        r = rule(WeightedDomain(Interval(0, 1), k); degree = 16, T)
        @test eltype(r) == T && eltype(weights(r)) == (k.kind === :exp ? Complex{T} : T)
        @test passed(verify(r))
        @test abs(integrate(x -> exp(β * x), r) - osc_exact(k.kind, β, 30, 0, 1)) < 50eps(T)
    end
    @test measure(WeightedDomain(Interval(0, 1), Oscillatory(2))) ≈ (exp(2im) - 1) / 2im
    @test measure(WeightedDomain(Interval(0, 1), Oscillatory(2, sin))) ≈ (1 - cos(2)) / 2
end

@testset "Filon rules: behaviour in ω and N" begin
    β = 7 // 10
    # the error falls like ω⁻² at a fixed number of points
    errs = map((100, 1000, 10^4)) do ω
        r = rule(WeightedDomain(Interval(0, 1), Oscillatory(ω)); degree = 8, digits = 30)
        abs(integral(r, β) - osc_exact(:exp, β, ω, 0, 1))
    end
    @test errs[2] < errs[1] / 50 && errs[3] < errs[2] / 50
    # a real kernel with a parity, on a centred interval, gains a degree
    r = rule(WeightedDomain(Interval(-1, 1), Oscillatory(3, cos)); degree = 8, digits = 30)
    @test degree(r) == 9 && passed(verify(r))
    r = rule(WeightedDomain(Interval(-1, 1), Oscillatory(3, sin)); degree = 9, digits = 30)
    @test degree(r) == 10 && passed(verify(r))
    r = rule(WeightedDomain(Interval(0, 2), Oscillatory(3, cos)); degree = 8, digits = 30)
    @test degree(r) == 8
    # near ω = 0 the rule tends to Gauss–Lobatto, and misses degree N + 1 by less than the
    # delivered precision can show; the certificate carries the miss
    r = rule(WeightedDomain(Interval(0, 1), Oscillatory(1 // 10^6)); degree = 8, digits = 30)
    v = verify(r)
    @test passed(v) && v.sharp === nothing && occursin("not resolvable", v.sharp_note)
    @test 0 < certificate(r).next_error < 1e-25
    # degree 0 still has both end points
    r = rule(WeightedDomain(Interval(0, 1), Oscillatory(5)); degree = 0, digits = 20)
    @test npoints(r) == 2 && degree(r) == 1 && passed(verify(r))
    @test any(c -> c.key == "Filon1928", provenance(r).citations)
end

@testset "oscillatory weights: arguments and display" begin
    @test sprint(show, Oscillatory(2)) == "Oscillatory(2)"
    @test sprint(show, Oscillatory(1 // 2, sin)) == "Oscillatory(1//2, sin)"
    @test_throws ArgumentError Oscillatory(0)
    @test_throws ArgumentError Oscillatory(Inf, cos)
    @test_throws MethodError Oscillatory(1, tan)
    @test_throws NoRuleError rule(WeightedDomain(Interval(0, 1), Oscillatory(3)); degree = 4, T = Rational{BigInt})
    @test [c.family for c in available(WeightedDomain(Interval(0, 1), Oscillatory(3)))] == ["Filon"]
    s = sprint(show, MIME"text/plain"(), rule(WeightedDomain(Interval(0, 1), Oscillatory(3)); degree = 4))
    @test occursin("weights   : complex", s)
end
