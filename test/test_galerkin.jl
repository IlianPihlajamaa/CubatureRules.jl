using CubatureRules, SauterSchwabQuadrature, StaticArrays, LinearAlgebra, Test

# Galerkin pairs (notes/singular-bem.md, stage 4) are left to SauterSchwabQuadrature.jl, which
# takes its one-dimensional rule as an argument and computes in that rule's arithmetic. Fed this
# package's Gauss–Legendre rules on [0, 1] at 60 digits, it integrates 1/|x − y| over pairs of
# triangles beyond Float64 accuracy. These tests guard that interplay and the vertex orderings
# its regularising maps expect: a common vertex first in both triangles, a common edge at
# positions 1 and 3 of both, in the same order.

gauss01(n) = (g = rule(GaussLegendre(), Interval(0, 1); degree = 2n - 1, digits = 60);
              collect(zip(first.(nodes(g)), weights(g))))

"`∫_{T1} ∫_{T2} 1/|x − y| dy dx` by SauterSchwabQuadrature, triangles parametrised as it expects."
function laplace_pair(t1, t2, strategy, n)
    J(t) = abs((t[1] - t[3])[1] * (t[2] - t[3])[2] - (t[1] - t[3])[2] * (t[2] - t[3])[1])
    x(t, u) = t[3] + u[1] * (t[1] - t[3]) + u[2] * (t[2] - t[3])
    return sauterschwab_parameterized((u, v) -> J(t1) * J(t2) / norm(x(t1, u) - x(t2, v)), strategy(gauss01(n)))
end

@testset "Galerkin pairs through SauterSchwabQuadrature.jl, at BigFloat" begin
    setprecision(BigFloat, 256) do
        P(x, y) = SVector{2,BigFloat}(x, y)
        # the same triangle twice, against the closed form for sides a, b, c and area A:
        # (4A²/3) Σ_cyc (1/a) ln(((a + b)² − c²)/(b² − (c − a)²))
        t = [P(2, 0), P(big(3) / 10, big(9) / 10), P(0, 0)]
        a, b, c = norm(t[2] - t[3]), norm(t[3] - t[1]), norm(t[1] - t[2])
        A = abs((t[1] - t[3])[1] * (t[2] - t[3])[2] - (t[1] - t[3])[2] * (t[2] - t[3])[1]) / 2
        term(a, b, c) = log(((a + b)^2 - c^2) / (b^2 - (c - a)^2)) / a
        exact = 4A^2 / 3 * (term(a, b, c) + term(b, c, a) + term(c, a, b))
        val = laplace_pair(t, t, CommonFace, 20)
        @test val isa BigFloat
        @test abs(val - exact) < 1e-16 * exact                # past Float64
        # a common vertex and a common edge: geometric convergence in the expected orderings
        tv1, tv2 = [P(0, 0), P(1, 0), P(0, 1)], [P(0, 0), P(-1, 0), P(0, -1)]
        @test abs(laplace_pair(tv1, tv2, CommonVertex, 16) - laplace_pair(tv1, tv2, CommonVertex, 12)) < 1e-12
        te1, te2 = [P(1, 0), P(0, 0), P(0, 1)], [P(1, 0), P(1, 1), P(0, 1)]
        @test abs(laplace_pair(te1, te2, CommonEdge, 16) - laplace_pair(te1, te2, CommonEdge, 12)) < 1e-11
        # well separated: a few points suffice
        tp = [P(2, 0), P(3, 0), P(2, 1)]
        @test abs(laplace_pair(tv1, tp, PositiveDistance, 12) - laplace_pair(tv1, tp, PositiveDistance, 8)) < 1e-15
    end
end
