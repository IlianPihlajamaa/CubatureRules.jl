using CubatureRules, Test
const CR = CubatureRules

# The scripts in scripts/ and benchmark/ are not run by the test suite, so a name that stops
# being exported breaks them silently. Twice already: verify_tables.jl called `degree_range`
# and construction.jl called `benchmark_construction` after both became public rather than
# exported, and the construction benchmark then failed on every CI run until noticed. This
# parses each script and flags calls to public, non-exported names that it neither imports
# nor defines itself. It does not run them.

"Names a script calls (functions and macros), imports from CubatureRules, and defines."
function script_names(src::AbstractString)
    called, imported, defined = Set{Symbol}(), Set{Symbol}(), Set{Symbol}()
    walk(x) = nothing
    function walk(x::Expr)
        if x.head in (:import, :using)
            for a in x.args
                a isa Expr && a.head === :(:) || continue
                a.args[1] == Expr(:., :CubatureRules) || continue
                for b in a.args[2:end]
                    b isa Expr && b.head === :. && push!(imported, b.args[end])
                    b isa Expr && b.head === :as && push!(imported, b.args[2])
                end
            end
        elseif x.head === :call && x.args[1] isa Symbol
            push!(called, x.args[1])
        elseif x.head === :macrocall && x.args[1] isa Symbol
            push!(called, Symbol(string(x.args[1])[2:end]))
        elseif x.head in (:function, :(=)) && x.args[1] isa Expr && x.args[1].head === :call &&
               x.args[1].args[1] isa Symbol
            push!(defined, x.args[1].args[1])
        elseif x.head === :(=) && x.args[1] isa Symbol
            push!(defined, x.args[1])
        end
        foreach(walk, x.args)
    end
    walk(Meta.parseall(src))
    return called, imported, defined
end

const PUBLIC_ONLY = Set(n for n in names(CR) if !Base.isexported(CR, n))

@testset "scripts use only names they can reach" begin
    root = pkgdir(CR)
    for dir in ("scripts", "benchmark"), f in sort(readdir(joinpath(root, dir); join = true))
        endswith(f, ".jl") || continue
        called, imported, defined = script_names(read(f, String))
        unreachable = sort!(collect(setdiff(intersect(called, PUBLIC_ONLY), imported, defined)))
        @testset "$(relpath(f, root))" begin
            @test isempty(unreachable)
        end
    end
    # and the guard itself catches the case it was written for
    called, imported, defined = script_names("using CubatureRules\nbenchmark_construction(Simplex{2}())")
    @test :benchmark_construction in setdiff(intersect(called, PUBLIC_ONLY), imported, defined)
    called, imported, defined = script_names("import CubatureRules: benchmark_construction\nbenchmark_construction(x)")
    @test isempty(setdiff(intersect(called, PUBLIC_ONLY), imported, defined))
end
