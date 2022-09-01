# Traces.jl : Graph canonical labeling and automorphism group computation

Simple wrapper for [`traces`](https://pallini.di.uniroma1.it/) of version 27r4  in Julia. 

## Installation
Require `gcc` and a POSIX style build environment. 

Add the package:
```julia
pkg.add("https://github.com/lidingxu/Traces.jl.git")
```


## Example usage


Convert a `Graph` g to a nauty's `SparseGraph` sp_g, and get its canonical labelling (permutation), automophism (generators), orbits:

```julia
using Traces
using Graphs
g = Graph()
# add vertices and edges to g ...
canon_label, automorphism, orbit_class = Traces.traces(g, true, true)
```


## API

Data structures
* `DEFAULTOPTIONS_TRACES` :  default options
* `TracesOptions` : defaul constructor of options
* `TracesStats` : defaul constructor of stats
* `SparseGraph`: defaul constructor of traces' `SparseGraph`
* `tracesreturn`: a struct containing return of traces, i.e., canonocial graph, generators, labels, partition, orbits, stats 

> **Note**:  the index of Julia and `Graph` start at 1, the index of C and `SparseGraph` start at 0. See the [user guide](https://pallini.di.uniroma1.it/Guide.html) of nauty and traces for details about their data structure and function call.


Interfaces
* `backend_traces`: a wrapper function call to C function `Traces` via ccall
* `backend_traces_with_automs`: a wrapper function call to C function `Traces_With_Automs` via ccall

Helper methods
* `to_sparse`: return a `SparseGraph` of `Graphs`, `SimpleGraphs` and `MetaGraphs`
* `traces`: return readeable orbits, canonocial graph, generators (in cyclic representation), of `Graphs`, `SimpleGraphs` and `MetaGraphs` give labels, parition

> **Warning**:  reuse of `SparseGraph` return by `to_sparse` may lead to a memory leak, because the memory allocated for the internal data structure may be released.