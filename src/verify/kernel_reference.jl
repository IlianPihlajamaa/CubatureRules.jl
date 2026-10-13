# Reference integrals for the kernels singular or near-singular at a point x₀ of a triangle,
# a curved triangle or a tetrahedron: the orthonormal polynomials integrated with the kernel
# in polar coordinates about x₀ (or about the point nearest it), independently of the
# constructions, which work in Duffy coordinates.

"""
    refined_integrals(integrals, levels, S, what) -> Vector{S}

`integrals(level)` at the successive `levels` until two agree to `precision(S) − 16` bits
(relative to the largest integral, or 1); the later of the two is returned. Throws, naming
`what`, if no two successive levels agree.
"""
function refined_integrals(integrals, levels, ::Type{S}, what) where {S}
    tol = ldexp(one(S), -(precision(S) - 16))
    prev = nothing
    for level in levels
        cur = integrals(level)
        prev !== nothing && maximum(abs, cur - prev) <= tol * max(one(S), maximum(abs, cur)) && return cur
        prev = cur
    end
    throw(ArgumentError("$what did not converge"))
end

# A triangle with the kernel 1/|y − x₀| is verified against the orthonormal Dubiner
# polynomials of the triangle, with their weighted integrals computed in polar coordinates
# about x₀ (see `patch_frame`): on each sub-triangle (x₀, p, q) the kernel cancels the polar
# Jacobian, the radial integral of a polynomial is exact with Gauss–Legendre, and the angle is
# taken along the edge pq, as τ = h sinh v, where the integrand is entire. The construction
# shares none of this: it works in the Duffy coordinates, with a Gauss rule of 1/√q(t) from a
# recurrence.
struct KernelTriangleBasis{S,D,L} <: VerificationBasis
    dub::DubinerBasis{S}
    A::SMatrix{2,D,S,L}             # reference coordinates ξ = A (y − v₀); a left inverse in space
    v0::SVector{D,S}
    integrals::Vector{S}
end
function evaluate!(B::KernelTriangleBasis{S,D}, y) where {S,D}
    ξ = B.A * (SVector{D,S}(ntuple(i -> S(y[i]), D)) - B.v0)
    dubiner!(B.dub.φ, B.dub.gx, B.dub.gy, B.dub.ws, ξ[1], ξ[2])
    # the gradient in y, for the tolerance on rounded nodes: Aᵀ ∇_ξ
    g = [sum(abs(B.A[1, i] * gx + B.A[2, i] * gy) for i in 1:D) for (gx, gy) in zip(B.dub.gx, B.dub.gy)]
    return B.dub.φ, g
end
blocks(B::KernelTriangleBasis) = blocks(B.dub)
exact_integrals(B::KernelTriangleBasis) = B.integrals

"""
    polar_kernel_integrals(dom, L, S, φ!; degree) -> Vector{S}

`∫_T φᵢ(y)/|y − x₀| dy` for `L` functions, `φ!(out, y)` writing their values at `y`, which are
polynomials of degree at most `degree`; in polar coordinates about the apex `c` of
[`kernel_geometry`](@ref), in arithmetic `S`. On each sub-triangle `(c, p, q)` the angle is
taken along the edge `pq`, `τ = h sinh v` from the foot of the perpendicular at distance `h`,
and the ray to the edge point `e` has length `R`; with `β = (e − c)·(c − x₀)/R` and
`H = |x₀ − c|` the sub-triangle contributes

    ∫ dv ∫₀¹ φ(c + ρ (e − c)) h ρ R / √(R² ρ² + 2 β R ρ + H²) dρ.

For `x₀` on the triangle (`H = 0`) the radial factor is `h`, and Gauss–Legendre in `ρ` is
exact; off it, `ρ = ρ* + κ sinh z` about the point of the ray nearest `x₀` turns the radial
factor into `h ρ dz`, entire in `z`. The angular integrand is entire in `v`. Both are refined
until two resolutions agree to the precision of `S`.

With `x₀` off the triangle the other point kernels are ordinary integrals too: the radial
factor `h ρ dz` is then multiplied by the kernel times `r`, `1/r²` for `1/r³` and
`(y − x₀)·e/r²` for the gradient kernel, with `r = R κ cosh z` along the ray.

How independent this is of the constructions: the radial direction is (they use Gauss rules
of the radial weight, from its moments), and so are the coordinates, polar about `c` against
Duffy's. The angular substitution is the one `DuffySinh` makes in its angular direction, set
up separately, along the edge, at its own resolution.
"""
function polar_kernel_integrals(dom::PointKernelDomain, L::Integer, ::Type{S}, φ!; degree::Integer = 0) where {S}
    g = kernel_geometry(dom)
    D = length(g.x0)
    x0, c = SVector{D,S}(g.x0), SVector{D,S}(g.c)
    H2 = S(g.H2)
    near = !iszero(g.H2)
    kernel = dom.weight              # not `w`: integrals() below assigns its Gauss weights to w
    kernel isa InverseDistance || near ||
        throw(ArgumentError("with x₀ on the triangle $(kernel) is a finite part: see polar_finitepart_integrals"))
    edir = kernel isa InverseDistanceGradient ? SVector{D,S}(map(S, kernel.e)) : nothing
    # the kernel times r, at the point y at distance r
    kr(y, r) = kernel isa InverseDistance ? one(S) : kernel isa InverseDistanceCubed ? 1 / r^2 : dot(y - x0, edir) / r^2
    nr = cld(degree + 1, 2) + 1                             # exact in ρ to degree 2nr − 1 > degree
    vals = zeros(S, L)
    function integrals(level)
        total = zeros(S, L)
        x, w = gauss_jacobi_work((2degree + 16) * 2^level, 0, 0, precision(S))[1:2]
        ρg, wg = gauss_jacobi_work(near ? (2degree + 16) * 2^level : nr, 0, 0, precision(S))[1:2]
        for (p, q) in g.patches
            f, u, h, τp, τq = patch_frame(c, p, q, S)
            # along the edge τ = τ* + κ sinh v: κ = h for x₀ on the triangle; off it τ* ± iκ are
            # the roots of |f + τ u − x₀|², which keeps the substitution to x₀'s distance
            ts = near ? -dot(f - x0, u) : zero(S)
            κa = near ? sqrt(max(dot(f - x0, f - x0) - ts^2, zero(S))) : h
            va, vb = asinh((τp - ts) / κa), asinh((τq - ts) / κa)
            vm, hv = (va + vb) / 2, (vb - va) / 2
            for (xk, wk) in zip(x, w)
                vv = vm + hv * xk
                e = f + (ts + κa * sinh(vv)) * u             # where the ray meets the edge
                if !near                                    # the radial factor is h: exact in ρ
                    for (ρj, wj) in zip(ρg, wg)
                        φ!(vals, c + (1 + ρj) / 2 * (e - c))
                        total .+= (hv * wk * h * wj / 2) .* vals
                    end
                    continue
                end
                # |c + ρ(e − c) − x₀|² = R²((ρ − ρ*)² + κ²): ρ = ρ* + κ sinh z turns the radial
                # factor h ρ R/√(…) dρ into h ρ dz, entire in z (ρ − ρ* = exp z when κ = 0);
                # dτ = R dv when κa = h, else κa cosh v dv
                R = norm(e - c)
                dτR = κa * cosh(vv) / R
                β = dot(e - c, c - x0) / R
                ρs = -β / R
                κ = sqrt(max(H2 - β^2, zero(S))) / R
                z0, z1 = iszero(κ) ? (log(-ρs), log(1 - ρs)) : (asinh(-ρs / κ), asinh((1 - ρs) / κ))
                zm, hz = (z0 + z1) / 2, (z1 - z0) / 2
                for (ρj, wj) in zip(ρg, wg)
                    z = zm + hz * ρj
                    ρ = iszero(κ) ? ρs + exp(z) : ρs + κ * sinh(z)
                    y = c + ρ * (e - c)
                    φ!(vals, y)
                    r = R * (iszero(κ) ? exp(z) : κ * cosh(z))
                    total .+= (hv * wk * dτR * hz * wj * h * ρ * kr(y, r)) .* vals
                end
            end
        end
        return total
    end
    return refined_integrals(integrals, 0:6, S, "the reference integrals for $(dom)")
end

"""
    polar_finitepart_integrals(dom, L, S, φ!; degree) -> Vector{S}

The finite parts with respect to `r = |y − x₀|` of `∫_T φᵢ(y) K(y) dy` for the kernel `K` of
a [`FinitePartDomain`](@ref), for `L` functions written by `φ!(out, y)` that are polynomials
of degree at most `degree`; in polar coordinates about `x₀`, in arithmetic `S`.

Along the ray to the edge point at distance `R`, `φ(x₀ + ρ ê) = Σ a_k (ρ/R)^k`, with the
`a_k` from the values at `degree + 1` points. The radial variable is `r` itself, so the
finite part along the ray is closed form: `FP ∫₀^R ρ^(k−2) dρ = R^(k−1)/(k − 1)`, and `ln R`
at `k = 1`, for `1/r³` (whose polar Jacobian leaves `ρ⁻²`); `FP ∫₀^R ρ^(k−1) dρ = R^k/k`, and
`ln R` at `k = 0`, times `ê·e` for the gradient kernel. The terms dropped (`1/ε`, `ln ε`) are
the ones the finite part drops. The angle is taken along each edge as `τ = h sinh v`, as in
[`polar_kernel_integrals`](@ref), and refined until two resolutions agree.

Independent of the construction, which works in the Duffy variable `s = r/√q(t)` and so
needs the `ln √q` corrections that the radial variable here does not.
"""
function polar_finitepart_integrals(dom::FinitePartDomain, L::Integer, ::Type{S}, φ!; degree::Integer = 0) where {S}
    g = kernel_geometry(dom)
    g.place === :near && throw(ArgumentError("the finite parts are defined here for x₀ on the triangle; $(dom) has it off"))
    D = length(g.x0)
    x0 = SVector{D,S}(g.x0)
    hyper = dom.weight isa InverseDistanceCubed
    e = hyper ? nothing : SVector{D,S}(map(S, dom.weight.e))
    K = Int(degree)
    σx = gauss_jacobi_work(K + 1, 0, 0, precision(S))[1]
    σ = (1 .+ S.(σx)) ./ 2
    Vinv = inv([σ[i]^k for i in 1:(K + 1), k in 0:K])
    vals = zeros(S, L)
    G = zeros(S, L, K + 1)
    function integrals(level)
        total = zeros(S, L)
        x, w = gauss_jacobi_work((2K + 16) * 2^level, 0, 0, precision(S))[1:2]
        for (p, q) in g.patches
            f, u, h, τp, τq = patch_frame(x0, p, q, S)
            va, vb = asinh(τp / h), asinh(τq / h)
            vm, hv = (va + vb) / 2, (vb - va) / 2
            for (xk, wk) in zip(x, w)
                v = vm + hv * xk
                edge = f + h * sinh(v) * u
                R = norm(edge - x0)
                ê = (edge - x0) / R
                for i in 1:(K + 1)
                    φ!(vals, x0 + σ[i] * R * ê)
                    copy_values!(view(G, :, i), vals)
                end
                a = G * transpose(Vinv)                      # a[l, k + 1]: coefficient of (ρ/R)^k
                dθ = hv * wk / cosh(v)                       # dθ = dv/cosh v along the edge
                c = hyper ? one(S) : dot(ê, e)
                for l in 1:L
                    acc = if hyper
                        -a[l, 1] + (K >= 1 ? a[l, 2] * log(R) : zero(S)) +
                        sum((a[l, k + 1] / (k - 1) for k in 2:K); init = zero(S))
                    else
                        a[l, 1] * log(R) + sum((a[l, k + 1] / k for k in 1:K); init = zero(S))
                    end
                    total[l] += dθ * c * (hyper ? acc / R : acc)
                end
            end
        end
        return total
    end
    return refined_integrals(integrals, 0:6, S, "the reference finite parts for $(dom)")
end

"""
    polar_curved_integrals(dom, L, S, φ!; degree) -> Vector{S}

`∫_T̂ φᵢ(ξ) J(ξ)/|χ(ξ) − χ(ξ₀)| dξ` for the [`CurvedInverseDistance`](@ref) weight of `dom`,
for `L` functions written by `φ!(out, ξ)`, in polar coordinates about `ξ₀` in arithmetic `S`:
`ρ dρ dθ · J/|χ(ξ₀ + ρ ê) − χ(ξ₀)| = J dρ dθ/|Dχ(ξ₀) ê + ρ Q(ê)|`, smooth in `ρ` (Gauss–Legendre
on each ray) and in the angle (taken along each edge as `τ = h sinh v`, as in
[`polar_kernel_integrals`](@ref)), both refined until two resolutions agree. Independent of the
construction, which works in Duffy coordinates with a Gauss rule of each ray's weight.

For `x₀` off the surface, about the point `ξc` nearest it: `ρ J dρ dθ/|D₀ + ρ Dχ ê + ρ² Q(ê)|`,
`D₀ = χ(ξc) − x₀`, with `ρ = ρ* + κ sinh z` from the quadratic part along each ray, as in
`polar_kernel_integrals`.
"""
function polar_curved_integrals(dom::CurvedKernelDomain, L::Integer, ::Type{S}, φ!; degree::Integer = 0) where {S}
    g = curved_geometry(dom; bits = precision(S) + 64)
    near = g.place === :near
    el = dom.weight.element
    X0 = SVector{2,S}(map(S, g.x0))
    A = jacobian_matrix(el, X0)
    D0 = near ? el(X0) - S.(g.X) : nothing
    vals = zeros(S, L)
    # the radial integral along one ray of length R, direction ê, with n Gauss–Legendre points
    function ray!(total, ê, R, x, w, scale)
        Ae, Qe = A * ê, second_term(el, ê)
        if !near
            for (xj, wj) in zip(x, w)
                ρ = R * (1 + xj) / 2
                ξ = X0 + ρ * ê
                φ!(vals, ξ)
                total .+= (scale * R / 2 * wj * area_element(el, ξ) / norm(Ae + ρ * Qe)) .* vals
            end
            return
        end
        q = dot(Ae, Ae)
        ρs = -dot(D0, Ae) / q
        κ2 = dot(D0, D0) / q - ρs^2
        κ = sqrt(max(κ2, zero(S)))
        z0, z1 = κ2 > 0 ? (asinh(-ρs / κ), asinh((R - ρs) / κ)) : (log(-ρs), log(R - ρs))
        zm, hz = (z0 + z1) / 2, (z1 - z0) / 2
        for (xj, wj) in zip(x, w)
            z = zm + hz * xj
            ρ, dρ = κ2 > 0 ? (ρs + κ * sinh(z), κ * cosh(z)) : (ρs + exp(z), exp(z))
            ξ = X0 + ρ * ê
            φ!(vals, ξ)
            total .+= (scale * hz * wj * dρ * ρ * area_element(el, ξ) / norm(D0 + ρ * Ae + ρ^2 * Qe)) .* vals
        end
    end
    function integrals(n)
        total = zeros(S, L)
        x, w = gauss_jacobi_work(n, 0, 0, precision(S))[1:2]
        for (p, q) in g.patches
            f, u, h, τp, τq = patch_frame(X0, p, q, S)
            # along the edge τ = τ* + κ sinh v, dθ = h dτ/(h² + τ²): on the surface κ = h; off
            # it τ* ± iκ are the roots of the linearised |χ(e) − x₀|², x₀'s distance included
            ts, κ = zero(S), h
            if near
                Dw, Au = D0 + A * (f - X0), A * u
                ts = -dot(Dw, Au) / dot(Au, Au)
                κ = sqrt(max(dot(Dw, Dw) / dot(Au, Au) - ts^2, zero(S)))
            end
            va, vb = asinh((τp - ts) / κ), asinh((τq - ts) / κ)
            vm, hv = (va + vb) / 2, (vb - va) / 2
            for (xk, wk) in zip(x, w)
                vv = vm + hv * xk
                τ = ts + κ * sinh(vv)
                e = f + τ * u
                R = norm(e - X0)
                ray!(total, (e - X0) / R, R, x, w, hv * wk * h * κ * cosh(vv) / (h^2 + τ^2))
            end
        end
        return total
    end
    ns = [2degree + 24]                                 # doubled while at most 2048
    while last(ns) <= 2048
        push!(ns, 2last(ns))
    end
    return refined_integrals(integrals, ns, S, "the reference integrals for $(dom)")
end

function verification_basis(dom::CurvedKernelDomain, n, ::Type{S}, exact) where {S}
    exact && throw(ArgumentError("rules for the curved kernel are not exact rationals"))
    K = maximum(_degrees(n))
    dub = DubinerBasis{S}(K; normalize = true)
    function φ!(out, ξ)
        dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, ξ[1], ξ[2])
        copy_values!(out, dub.φ)                  # by value: dubiner! writes into dub.φ in place
    end
    integrals = polar_curved_integrals(dom, dubiner_length(K), S, φ!; degree = K)
    A = SMatrix{2,2,S,4}(1, 0, 0, 1)
    return KernelTriangleBasis{S,2,4}(dub, A, SVector{2,S}(0, 0), integrals),
           "orthonormal Dubiner in ξ, against their integrals with the weight (polar coordinates about ξ₀, " *
           "Gauss–Legendre along each ray and each edge)"
end

# A tetrahedron with the kernel 1/|y − x₀| is verified against the orthonormal tetrahedral
# Dubiner polynomials, with their weighted integrals over each cone from x₀: Gauss–Legendre in
# the radial variable (exact, the integrand being a polynomial there) and, on the face, the
# polar quadrature of `polar_kernel_integrals` about the face's point nearest x₀. The
# construction shares the cones and nothing on the faces, where it uses DuffySinh.
struct KernelTetBasis{S,B<:SimplexBasis{3,S}} <: VerificationBasis
    basis::B
    A::SMatrix{3,3,S,9}             # reference coordinates ξ = A (y − v₀)
    v0::SVector{3,S}
    integrals::Vector{S}
end
function evaluate!(B::KernelTetBasis{S}, y) where {S}
    ξ = B.A * (SVector{3,S}(S(y[1]), S(y[2]), S(y[3])) - B.v0)
    φ, G = evaluate!(B.basis, ξ)
    # the gradient in y, for the tolerance on rounded nodes: Aᵀ ∇_ξ
    return φ, [sum(abs(B.A[1, i] * G[k, 1] + B.A[2, i] * G[k, 2] + B.A[3, i] * G[k, 3]) for i in 1:3) for k in axes(G, 1)]
end
blocks(B::KernelTetBasis) = degree_blocks(B.basis)
exact_integrals(B::KernelTetBasis) = B.integrals

function verification_basis(dom::InverseDistanceTetDomain, n, ::Type{S}, exact) where {S}
    exact && throw(ArgumentError("rules for the kernel 1/|y − x₀| are not exact rationals"))
    K = maximum(_degrees(n))
    v = [SVector{3,S}(map(S, p)) for p in vertices(dom.base)]
    A = inv(hcat(v[2] - v[1], v[3] - v[1], v[4] - v[1]))
    basis = SimplexBasis{3,S}(K)
    L = basis_length(basis)
    X0 = SVector{3,S}(map(S, tet_kernel_geometry(dom).x0))
    sx, sw = gauss_jacobi_work(cld(K + 2, 2), 0, 0, precision(S))[1:2]   # exact to degree K + 1 in s
    integrals = zeros(S, L)
    for (fd, h) in cone_faces(dom, S)
        function ψ!(out, z)
            fill!(out, zero(S))
            for (si, σi) in zip(sx, sw)
                s = (1 + S(si)) / 2
                φ, _ = evaluate!(basis, A * (X0 + s * (SVector{3,S}(z) - X0) - v[1]); gradient = false)
                out .+= (S(σi) / 2 * s) .* φ
            end
        end
        integrals .+= h .* polar_kernel_integrals(fd, L, S, ψ!; degree = K)
    end
    return KernelTetBasis{S,typeof(basis)}(basis, A, v[1], integrals),
           "orthonormal tetrahedral Dubiner, against their integrals with the kernel (cones from x₀ over the faces: " *
           "Gauss–Legendre along each cone, the face in polar coordinates about its point nearest x₀)"
end

function verification_basis(dom::FinitePartDomain, n, ::Type{S}, exact) where {S}
    exact && throw(ArgumentError("rules for the finite-part kernels are not exact rationals"))
    K = maximum(_degrees(n))
    D = length(dom.weight.x0)
    v = [SVector{D,S}(map(S, p)) for p in vertices(dom.base)]
    E = hcat(v[2] - v[1], v[3] - v[1])
    A = inv(E' * E) * E'
    dub = DubinerBasis{S}(K; normalize = true)
    function φ!(out, y)
        ξ = A * (y - v[1])
        dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, ξ[1], ξ[2])
        copy_values!(out, dub.φ)                  # by value: dubiner! writes into dub.φ in place
    end
    integrals = kernel_geometry(dom).place === :near ? polar_kernel_integrals(dom, dubiner_length(K), S, φ!; degree = K) :
                polar_finitepart_integrals(dom, dubiner_length(K), S, φ!; degree = K)
    return KernelTriangleBasis{S,D,2D}(dub, A, v[1], integrals),
           "orthonormal Dubiner, against their finite parts with the kernel (polar coordinates about x₀, " *
           "the finite part along each ray in closed form)"
end

function verification_basis(dom::InverseDistanceDomain, n, ::Type{S}, exact) where {S}
    exact && throw(ArgumentError("rules for the kernel 1/|y − x₀| are not exact rationals"))
    K = maximum(_degrees(n))
    D = length(dom.weight.x0)
    v = [SVector{D,S}(map(S, p)) for p in vertices(dom.base)]
    E = hcat(v[2] - v[1], v[3] - v[1])
    A = inv(E' * E) * E'                                    # E⁻¹ in the plane, its left inverse in space
    dub = DubinerBasis{S}(K; normalize = true)
    function φ!(out, y)
        ξ = A * (y - v[1])
        dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, ξ[1], ξ[2])
        copy_values!(out, dub.φ)                  # by value: dubiner! writes into dub.φ in place
    end
    integrals = polar_kernel_integrals(dom, dubiner_length(K), S, φ!; degree = K)
    return KernelTriangleBasis{S,D,2D}(dub, A, v[1], integrals),
           "orthonormal Dubiner, against their integrals with the kernel (polar coordinates about " *
           "the nearest point of the triangle, Gauss–Legendre along each edge and graded along each ray)"
end
