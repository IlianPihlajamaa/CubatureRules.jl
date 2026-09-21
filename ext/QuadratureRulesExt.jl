# Makes the QuadratureRules.jl-backed families constructible. Loading this package is all
# that is needed: the registry discovers families with `subtypes` at call time, so
# `available(Interval(); degree = …)` lists Lobatto, Radau and Clenshaw–Curtis from the
# moment QuadratureRules is loaded, with no registration step (PLAN §2.5).
module QuadratureRulesExt

using CubatureRules
using QuadratureRules
import CubatureRules: candidates, upstream_rule, Lobatto, Radau, ClenshawCurtis, Interval,
                      PolynomialDegree, isreference

const SYMMETRIC = QuadratureRules.SymmetricInterval()

candidates(::Type{Lobatto}, dom::Interval, ::PolynomialDegree) = isreference(dom) ? [Lobatto()] : Lobatto[]
candidates(::Type{Radau}, dom::Interval, ::PolynomialDegree) =
    isreference(dom) ? [Radau(:right), Radau(:left)] : Radau[]
candidates(::Type{ClenshawCurtis}, dom::Interval, ::PolynomialDegree) =
    isreference(dom) ? [ClenshawCurtis()] : ClenshawCurtis[]

upstream_rule(::Lobatto, n) = (QuadratureRules.lobatto_legendre_nodes(BigFloat, n; interval = SYMMETRIC),
                               QuadratureRules.lobatto_legendre_weights(BigFloat, n; interval = SYMMETRIC))
upstream_rule(f::Radau, n) = (QuadratureRules.radau_legendre_nodes(BigFloat, n, Val(f.side); interval = SYMMETRIC),
                              QuadratureRules.radau_legendre_weights(BigFloat, n, Val(f.side); interval = SYMMETRIC))
upstream_rule(::ClenshawCurtis, n) = (QuadratureRules.clenshaw_curtis_nodes(BigFloat, n; interval = SYMMETRIC),
                                      QuadratureRules.clenshaw_curtis_weights(BigFloat, n; interval = SYMMETRIC))

end
