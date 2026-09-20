# The v0.2 acceptance test (PLAN §3.2): integrating over a large mesh with a package rule
# must cost about what a hand-written FEM kernel costs.
#
#     julia --project -t auto benchmark/mesh.jl [ncells]
#
# The hand-written kernel has the same 6-point degree-4 rule inlined as literal constants and
# does the affine map by hand — what a FEM code would write if it did not use this package.

using CubatureRules, Printf, StaticArrays
using Base.Threads: nthreads

"A structured mesh of `2m²` triangles covering the unit square."
function square_mesh(m::Int)
    cells = Simplex{2,Float64,3}[]
    h = 1 / m
    for i in 0:(m - 1), j in 0:(m - 1)
        x, y = i * h, j * h
        push!(cells, Simplex((x, y), (x + h, y), (x + h, y + h)))
        push!(cells, Simplex((x, y), (x + h, y + h), (x, y + h)))
    end
    return cells
end

integrand(p) = exp(-p[1] * p[2]) * (1 + p[1]^2)

"The same rule, inlined by hand: barycentric coordinates and weights as literals."
function handwritten(f::F, cells) where {F}
    # degree-4 rule on the reference triangle, two 3-point orbits (Dunavant/Strang values)
    a1, w1 = 0.445948490915965, 0.111690794839005
    a2, w2 = 0.091576213509771, 0.054975871827661
    bary = ((a1, a1), (1 - 2a1, a1), (a1, 1 - 2a1), (a2, a2), (1 - 2a2, a2), (a2, 1 - 2a2))
    ws = (w1, w1, w1, w2, w2, w2)
    acc = 0.0
    @inbounds for c in cells
        v = c.vertices
        e1 = v[2] - v[1]
        e2 = v[3] - v[1]
        detJ = abs(e1[1] * e2[2] - e1[2] * e2[1])
        s = 0.0
        for k in 1:6
            ξ, η = bary[k]
            s += ws[k] * f(v[1] + ξ * e1 + η * e2)
        end
        acc += detJ * s
    end
    return acc
end

function main(ncells::Int)
    m = max(1, round(Int, sqrt(ncells / 2)))
    cells = square_mesh(m)
    r = static(rule(Simplex{2}(); degree = 4))
    @printf("%d cells, %d-point rule, %d threads\n", length(cells), npoints(r), nthreads())

    a = integrate(integrand, r, cells)                    # warm up
    b = handwritten(integrand, cells)
    c = integrate(integrand, r, cells; threaded = true)
    @printf("values agree: package %.15e, hand-written %.15e, threaded %.15e\n", a, b, c)

    t1 = @elapsed integrate(integrand, r, cells)
    t2 = @elapsed handwritten(integrand, cells)
    t3 = @elapsed integrate(integrand, r, cells; threaded = true)
    alloc = @allocated integrate(integrand, r, cells)
    @printf("package        : %.3f s  (%d bytes allocated)\n", t1, alloc)
    @printf("hand-written   : %.3f s\n", t2)
    @printf("package/hand   : %.2f×\n", t1 / t2)
    @printf("threaded       : %.3f s  (%.2f× the serial package path)\n", t3, t1 / t3)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? 1_000_000 : parse(Int, ARGS[1]))
end
