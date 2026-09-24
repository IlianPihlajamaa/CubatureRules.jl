# Display (PLAN §3.5): make correctness inspectable rather than a matter of trust. The
# residual line always names which check it reports — the certificate's defining-equation
# residual, never to be confused with verification.

function Base.show(io::IO, r::QuadratureRule{D,T}) where {D,T}
    print(io, "QuadratureRule{", D, ",", T, "}(", family(r), ", ", npoints(r), " points, ",
          describe(r.exactness), ", on ", r.domain, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::QuadratureRule{D,T}) where {D,T}
    p = r.provenance
    c = r.certificate
    println(io, "QuadratureRule{", D, ",", T, "} on ", r.domain)
    how = p.derivation isa Derived ? "derived" : "seeded"
    refined = c !== nothing && c.iterations > 0 ? ", Newton-refined" : ""
    println(io, "  family    : ", p.family, " (", how, refined, ")")
    claimtxt = r.exactness isa PolynomialDegree ? "polynomial degree $(r.exactness.d) (claimed; `check(rule)` verifies)" :
               describe(r.exactness)
    println(io, "  exactness : ", claimtxt)
    inside = all(x -> isinterior(x, r.domain), r.nodes)
    println(io, "  points    : ", npoints(r), inside ? ", all interior" : ", not all interior")
    pos = all(>(0), r.weights)
    println(io, "  weights   : ", pos ? "all positive" : "some negative", ", Σw = ", _short(sum(r.weights)))
    if c !== nothing
        prec = _is_exact_type(T) ? "exact rational" : "$(c.digits) digits"
        extra = String[]
        c.guard_digits > 0 && push!(extra, "guard $(c.guard_digits)")
        c.cond != 1 && push!(extra, @sprintf("est. cond(J) = %.1e", c.cond))
        println(io, "  precision : ", prec, isempty(extra) ? "" : " (" * join(extra, ", ") * ")")
        rdig = c.residual_bits > 0 ? " at $(floor(Int, c.residual_bits * log10(2)))-digit arithmetic" : ""
        println(io, "  residual  : ", @sprintf("%.1e", c.residual), "  (defining equations", rdig, ": ", c.equations, ")")
    end
    if !isempty(p.citations)
        cit = first(p.citations)
        print(io, "  reference : ", _short_author(cit), " (", cit.year, ")", isempty(cit.doi) ? "" : ", doi:" * cit.doi)
    else
        print(io, "  reference : —")
    end
end

_short(x::Rational) = string(x)
_short(x::BigFloat) = string(BigFloat(x; precision = 100))[1:min(end, 32)]
_short(x) = string(x)

function _short_author(c::Citation)
    last(n) = split(n)[end]
    length(c.authors) == 1 && return last(c.authors[1])
    length(c.authors) == 2 && return last(c.authors[1]) * " & " * last(c.authors[2])
    return last(c.authors[1]) * " et al."
end

function Base.show(io::IO, ::MIME"text/plain", c::Certificate)
    println(io, "Certificate (defining equations — not a verification)")
    println(io, "  equations  : ", c.equations)
    println(io, @sprintf("  residual   : %.2e at %d bits", c.residual, c.residual_bits))
    println(io, "  digits     : ", c.digits == typemax(Int) ? "exact" : c.digits, "  (guard ", c.guard_digits, ")")
    print(io, @sprintf("  cond(J)    : %.2e,  iterations: %d", c.cond, c.iterations))
    c.next_error === nothing || print(io, @sprintf("\n  next degree: missed by %.2e at working precision", c.next_error))
end

function Base.show(io::IO, ::MIME"text/plain", p::Provenance)
    println(io, "Provenance: ", p.family, " (", p.derivation isa Derived ? "derived" : "seeded", ")")
    for s in p.path
        println(io, "  • ", s)
    end
    println(io, "  seed      : ", p.seed_source)
    println(io, "  licence   : ", p.license)
    isempty(p.selection) || println(io, "  selection : ", p.selection)
    print(io, "  cites     : ", join((c.key for c in p.citations), ", "))
end
