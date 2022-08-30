# Traces.jl : Graph canonical labeling and automorphism group computation

Simple wrapper for [`traces`](https://pallini.di.uniroma1.it/) of version 27r4  in Julia. 

## Installation
Require `gcc` and a POSIX style build environment. 

Add the package:
```julia
pkg.add("https://github.com/lidingxu/Traces.jl.git")
```


## Example usage


Convert a `Graph` g to a nauty's `SparseGraph` sp_g:

```julia
using Traces
using Graphs
g = Graph()
# add vertices and edges to g ...
sp_g = to_sparse(g)
```


Compute orbits of sp_g (g):

```julia
orbits_g = orbits(sp_g)
```


## API

Data structures
* `DEFAULTOPTIONS_TRACES` :  default options
* `TracesOptions` : defaul constructor of options
* `TracesStats` : defaul constructor of stats
* `SparseGraph`: defaul constructor of traces' `SparseGraph`
* `tracesreturn`: a struct containing return of traces, i.e., canonocial graph, labels, partition, orbits, stats 

> **Note**:  the index of Julia and `Graph` start at 1, the index of C and `SparseGraph` start at 0. See the [user guide](https://pallini.di.uniroma1.it/Guide.html) of nauty and traces for details about their data structure and function call.
 **Warning**:  reuse of `SparseGraph` return may lead to memory leak, because the memory allocated by this data structure may be released.

Interfaces
* `traces`: a wrapper function call to C function traces via ccall

Helper methods
* `to_sparse`: return a `SparseGraph` of `Graphs`, `SimpleGraphs` and `MetaGraphs`
* `orbits`: return a dictionary of sets (orbits) of a `SparseGraph` 
* `to_label`: return array-like labelling, parition of a dictionary of sets (labels)

