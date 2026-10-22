# Record `next_error` in the sphere table: how far each rule misses degree d + 1, measured by
# `check` on the rule refined to 90 digits.
#
#     julia --project scripts/record_next_error.jl [degree ...]
#
# Without degrees, every entry from degree 19 on. In Float64 these rules come too close to
# the next degree for `check` to resolve the miss (at degree 131 it is 1e-39); recorded in the
# table, it lets a shipped rule report sharpness as not resolvable rather than fail (see
# `next_error` in `Certificate`). The degrees 3-17 miss by far more than Float64 resolves.
#
# A rule that is not exact and sharp at 90 digits stops the script. The table is rewritten
# line by line: the `next_error` line of each entry is replaced, or inserted after
# `min_weight` (the keys are sorted), and every other line is left as it is.

using CubatureRules, Printf
const CR = CubatureRules

const FILE = joinpath(@__DIR__, "..", "src", "data", "lebedev_seeds.toml")

function next_error(d)
    r = rule(LebedevRule(), Sphere{3}(); degree = d, digits = 90)
    v = check(r)
    (v.exact && v.sharp === true) ||
        error("degree $d: not exact and sharp at 90 digits (exact = $(v.exact), sharp = $(v.sharp))")
    return Float64(v.sharp_residual)
end

function record!(lines, d, value)
    i = findfirst(==("degree = $d"), lines)
    i === nothing && error("no entry of degree $d in $FILE")
    j = findnext(==("[[rule]]"), lines, i)
    block = i:(j === nothing ? length(lines) : j - 1)
    k = findfirst(l -> startswith(l, "next_error = "), lines[block])
    if k === nothing
        m = findfirst(l -> startswith(l, "min_weight = "), lines[block])
        insert!(lines, first(block) + m, "next_error = $(repr(value))")
    else
        lines[first(block) + k - 1] = "next_error = $(repr(value))"
    end
    return lines
end

if abspath(PROGRAM_FILE) == @__FILE__
    degrees = isempty(ARGS) ? [e.degree for e in CR.lebedev_entries() if e.degree >= 19] : parse.(Int, ARGS)
    lines = readlines(FILE)
    for d in degrees
        t = @elapsed miss = next_error(d)
        @printf("degree %3d: misses degree %d by %.3e  [%.1f s]\n", d, d + 1, miss, t)
        flush(stdout)
        record!(lines, d, miss)
        write(FILE, join(lines, "\n") * "\n")      # after every degree, so a long run can stop
    end
end
