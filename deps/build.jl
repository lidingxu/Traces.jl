@info "Building traces..."
using Libdl

function packagedir(pkg::AbstractString)
    if endswith(pkg, ".jl")
        pkg = pkg[1:end-3]
    end
    pkgdir = Base.find_package(pkg)
    isnothing(pkgdir) && throw(ErrorException("package \"$pkg\" not found"))              
    pkgdir = abspath(joinpath(dirname(pkgdir), ".."))
    return pkgdir
end

depsdir = joinpath(packagedir("Traces"), "deps") 
tracesdir = joinpath(depsdir, "nauty27r4")
tracesfiles = joinpath.(tracesdir, ["traces", "nauty","nautil", "naugraph", "schreier","naurng", "gtools"] .* ".c")
traceswrapper = joinpath(depsdir, "mintraceswrap.c") 
traceslib = joinpath(depsdir, "mintraceswrap." * Libdl.dlext) 

run(`gcc -DWORDSIZE=64 -DMAXN=0 -O4 -o $traceslib $traceswrapper $tracesfiles -shared -fPIC -I $tracesdir`)
