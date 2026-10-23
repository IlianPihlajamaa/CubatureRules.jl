# Construction-cost regression benchmark, run in CI on every commit.
#     julia --project benchmark/construction.jl [out.csv [previous.csv]]
# With the previous run's results, a case more than 3× slower or allocating more than 2×
# what it did fails the run; tiny cases, where the runner's noise dominates, are not judged.
using CubatureRules, Printf
import CubatureRules: benchmark_construction     # public, not exported

out = length(ARGS) >= 1 ? ARGS[1] : nothing
previous = length(ARGS) >= 2 && isfile(ARGS[2]) ? ARGS[2] : nothing
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

if previous !== nothing
    # case,family,degree → (seconds, bytes) of the previous run
    before = Dict{Tuple{String,String,Int},Tuple{Float64,Float64}}()
    for line in Iterators.drop(eachline(previous), 1)
        m = match(r"^\"(.*)\",(.*),(\d+),(\d+),([^,]+),([^,]+)$", line)
        m === nothing && continue
        before[(m[1], m[2], parse(Int, m[3]))] = (parse(Float64, m[5]), parse(Float64, m[6]))
    end
    worse = String[]
    for (label, r) in rows
        key = (label, r.family, r.degree)
        haskey(before, key) || continue
        s0, b0 = before[key]
        slow = r.seconds > 0.05 && r.seconds > 3s0
        fat = r.bytes > 1e6 && r.bytes > 2b0
        (slow || fat) && push!(worse, @sprintf("%s, %s, degree %d: %.3f s (was %.3f), %d bytes (was %d)",
                                               label, r.family, r.degree, r.seconds, s0, r.bytes, Int(b0)))
    end
    println("\ncompared with the previous run: ", length(before), " cases, ", length(worse), " regressed")
    foreach(w -> println("  ", w), worse)
    isempty(worse) || error("construction got slower or allocates more; see above")
end
