# A precompile workload over the common paths, so that a first `rule(...)` in a new session
# does not spend its time compiling. Before this, the first `Float64` triangle rule took 22 s,
# 98% of it compilation, for a rule that is then returned from a table in milliseconds.
#
# The workload runs while the package is precompiled. The caches it fills are safe to keep:
# the invariant bases and invariant moments are deterministic, and the family list is
# emptied in `__init__`.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    domains = (Simplex{2}(), Simplex{3}(), Orthotope{2}(), Orthotope{3}(), Sphere{3}(), Interval(), Ball{3}())
    @compile_workload begin
        for dom in domains
            r = rule(dom; degree = 5)
            integrate(x -> sum(x), r)
            show(IOBuffer(), MIME"text/plain"(), r)
        end
        r = rule(Simplex{2}(); degree = 4)
        check(r)
        integrate(x -> x[1] * x[2], r, [Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))])
        # one refinement above Float64 per kind of family: seeded, and Gauss on an interval
        rule(Simplex{2}(); degree = 4, digits = 30)
        rule(Interval(); degree = 7, digits = 30)
        available(Simplex{2}())
    end
end
