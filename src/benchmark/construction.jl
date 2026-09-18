# Construction-only benchmarking (PLAN §4.6): the package's own performance regression
# test, tracked from v0.1.

"""
    benchmark_construction(domain; degrees = 1:2:21, T = Float64, digits = nothing,
                           families = :all, warmup = true)

Time rule construction for every applicable family (or the given vector of family
instances) at each degree, with no integrand. Each family is warmed up at its lowest
degree first, so compilation is not timed; each point is a single sample, recorded as such.

Returns a vector of NamedTuples `(family, degree, npoints, seconds, bytes, status)`, where
`status` is `:completed` or `:unsupported`.
"""
function benchmark_construction(dom::Domain; degrees = 1:2:21, T = nothing, digits = nothing,
                                families = :all, warmup::Bool = true)
    Tout, _ = resolve_precision(T, digits)
    insts = families === :all ? applicable_instances(dom) : collect(families)
    rows = NamedTuple[]
    ref = reference(dom)
    for f in insts
        rg = degree_range(f, ref)
        if !supports_type(f, Tout)
            push!(rows, (family = describe_family(f), degree = -1, npoints = 0, seconds = NaN, bytes = 0,
                         status = :unsupported, samples = 0))
            continue
        end
        warmup && !isempty(rg) && rule(f, dom; degree = first(rg), T, digits)
        for d in degrees
            if !(d in rg)
                push!(rows, (family = describe_family(f), degree = d, npoints = 0, seconds = NaN, bytes = 0,
                             status = :unsupported, samples = 0))
                continue
            end
            stats = @timed rule(f, dom; degree = d, T, digits)
            push!(rows, (family = describe_family(f), degree = d, npoints = npoints(stats.value),
                         seconds = stats.time, bytes = stats.bytes, status = :completed, samples = 1))
        end
    end
    return rows
end
