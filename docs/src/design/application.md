# Application and transport

## Integration

`integrate(f, r)` computes `Σᵢ wᵢ f(xᵢ)`. In one dimension `f` receives a number, otherwise
an `SVector` of coordinates.

The sum starts from the first term `w₁ f(x₁)` rather than from a zero of some fixed type. The
result therefore has whatever type `f` returns: a complex number, a vector, a `StaticArray`,
or a quantity with units. Nothing is allocated if the nodes are `isbits` and `f` does not
allocate.

The integrand is typed `f::F` in every method, which makes Julia specialise on it even
though it is only passed on. Without this, the call would go through dynamic dispatch.

`integrate(f, r; batch = true)` calls `f` once with all nodes, as a `D × N` matrix (a vector
in one dimension), and expects the `N` values back. This suits integrands that are cheaper to
evaluate on many points at once.

## Mapping to another domain

`integrate(f, r, dom)` maps the nodes of a reference rule `r` onto the domain `dom` as they
are used, and multiplies by the Jacobian determinant. No new rule is built.

`integrate(f, r, cells)` does the same for each cell in a vector of domains and adds the
results, which is the usual way to integrate over a mesh. With `threaded = true` each thread
keeps its own accumulator, and the per-thread sums are added in thread order, so the result
depends on the number of threads but not on how the cells were scheduled.

[`map_to`](@ref)`(r, dom)` returns the mapped rule as a new `QuadratureRule`, with the step
added to its provenance.

Supported maps are the affine maps between simplices, intervals and boxes, and similarities
(translation and scaling) between spheres and balls. All of them preserve polynomial
degree.

## Subdivision

[`subdivide`](@ref)`(r, dom, n)` splits the domain into `n^D` congruent pieces (for simplices,
the standard subdivision into `n^D` smaller simplices), maps the rule onto each and
concatenates the results. The degree is preserved.

## Nonlinear transformations

[`transform`](@ref)`(r, φ, Jφ)` applies a change of variables `φ` with Jacobian `Jφ`. The
nodes become `φ(xᵢ)` and the weights `wᵢ |det Jφ(xᵢ)|`. Because `φ` is not affine, the
exactness claim is replaced by [`NoClaim`](@ref), unless the caller passes `claim`, in which
case the provenance records the claim as asserted by the caller.

[`duffy`](@ref)`(r; power = q)` is a specific transformation for simplices: it moves the
nodes towards the vertex at the origin by `x ↦ u^{q-1} x` with `u = Σxᵢ`, which cancels a
singularity of the form `r^{-(q-1)D}` at that vertex.

## Why there is only one rule type

Earlier versions also had a `static` form that stored nodes and weights as `SVector`s, with
the number of points in the type, and unrolled the sum at compile time. It was removed in
v0.5 after measurement on triangle rules of 7 to 445 points:

- with a nearly free integrand it was 1.0–1.4× faster in most cases, 1.9× only at 7 points
  (about 4 ns), and 0.84× (slower) at 445 points;
- with an integrand such as `exp` the difference was within noise;
- compiling each new size took from 0.2 s at 7 points to 2.5 s at 445 points;
- integration with the ordinary rule already allocates nothing.

Code that needs an `isbits` rule, such as a GPU kernel, can build one directly with
`SVector{N}(nodes(r))` and `SVector{N}(weights(r))`.

## Tensor products of rules

Two rules can be combined into a rule on the product domain with `r₁ ⊗ r₂`. The operator is
public but not exported, because other packages also define `⊗`; bring it in with
`using CubatureRules: ⊗`. For boxes, [`TensorProduct`](@ref) does the same through the
selector.
