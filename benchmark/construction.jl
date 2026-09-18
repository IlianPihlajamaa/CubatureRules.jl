# Construction-cost regression benchmark (PLAN §4.6), run in CI on every commit.
#     julia --project benchmark/construction.jl [out.csv]
using CubatureRules, Printf

out = isempty(ARGS) ? nothing : ARGS[1]
cases = [
    ("Simplex{2}, Float64", Simplex{2}(), (degrees = 1:2:19, T = Float64)),
    ("Simplex{2}, 100 digits", Simplex{2}(), (degrees = [5, 10, 15, 20], digits = 100)),
    ("Simplex{2}, 200 digits", Simplex{2}(), (degrees = [10, 20], digits = 200)),
    ("Simplex{3}, 50 digits", Simplex{3}(), (degrees = [3, 7, 11], digits = 50)),
]
rows = []
for (label, dom, kw) in cases
    for r in benchmark_construction(dom; kw...)
        r.status === :completed || continue
        push!(rows, (label, r))
        @printf("%-24s %-28s degree %2d  %5d points  %8.3f s  %10d bytes\n",
                label, r.family, r.degree, r.npoints, r.seconds, r.bytes)
    end
end
if out !== nothing
    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "case,family,degree,npoints,seconds,bytes")
        for (label, r) in rows
            println(io, "\"", label, "\",", r.family, ",", r.degree, ",", r.npoints, ",", r.seconds, ",", r.bytes)
        end
    end
end
