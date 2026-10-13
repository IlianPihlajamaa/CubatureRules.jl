# The symmetry groups a rule can claim, as types: a family states one in its `properties` and
# its `Provenance`, and verification checks invariance under it by dispatch
# (`check_invariance`, src/verify/symmetry.jl). Each acts on the reference domain's
# coordinates by permutations of one block of them, with or without sign changes.

"""
    SymmetryGroup

The symmetry group a rule claims to be invariant under (`provenance(r).symmetry`):
[`NoSymmetry`](@ref), [`Reflection`](@ref) (`x ↦ −x` on the interval), [`SymmetricGroup`](@ref)`{N}`
(the permutations of a simplex's `N` barycentric coordinates), [`Oh`](@ref) (the 48 signed
permutations of three coordinates: the cube, the sphere, the ball), [`D4`](@ref) (the 8 of
two: the square, the disk), [`B4`](@ref) (the 384 of four: the 4-cube), [`C4v`](@ref) (the
square's 8 on the pyramid's base coordinates) and [`D3h`](@ref) (the wedge's 12: the
triangle's 6 times `z ↦ −z`). Printed by their usual names.
"""
abstract type SymmetryGroup end

"No symmetry claimed."
struct NoSymmetry <: SymmetryGroup end
"`x ↦ −x` on the reference interval."
struct Reflection <: SymmetryGroup end
"The permutations of the `N` barycentric coordinates of a simplex (`S₃` on a triangle)."
struct SymmetricGroup{N} <: SymmetryGroup end
"The 48 signed permutations of three Cartesian coordinates (the octahedral group)."
struct Oh <: SymmetryGroup end
"The 8 signed permutations of two Cartesian coordinates (the square's group)."
struct D4 <: SymmetryGroup end
"The 384 signed permutations of four Cartesian coordinates (the hyperoctahedral group of the 4-cube)."
struct B4 <: SymmetryGroup end
"The square's 8 symmetries on the first two coordinates of the pyramid, the height fixed."
struct C4v <: SymmetryGroup end
"The wedge's 12: the permutations of the triangle's barycentric coordinates, times `z ↦ −z`."
struct D3h <: SymmetryGroup end

symmetry_name(::NoSymmetry) = "none"
symmetry_name(::Reflection) = "reflection"
symmetry_name(::SymmetricGroup{N}) where {N} = "S$N"
symmetry_name(::Oh) = "Oh"
symmetry_name(::D4) = "D4"
symmetry_name(::B4) = "B4"
symmetry_name(::C4v) = "C4v"
symmetry_name(::D3h) = "D3h"
Base.show(io::IO, g::SymmetryGroup) = print(io, symmetry_name(g))
