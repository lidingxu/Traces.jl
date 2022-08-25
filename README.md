# Traces.jl

Simple wrapper for using `traces` of version 27r4, a graph isomorphism package, with `SimpleGraphs`, `Graphs` and `MetaGraphs`,  in Julia. Requires `gcc` and a POSIX style build environment. 

## installation
```julia
pkg.add(")
```

## Example usage


Graphs, SimpleGraphs and MetaGraphs interface, convert a graph g to a nauty Sparsegraph sp_g:

```julia
using Traces
using Graphs
g = Graph()
# add vertices and edges to g ...
sp_g = to_sparse(g)
```

Compute orbits:

```julia
orbits(sp_g)
```

