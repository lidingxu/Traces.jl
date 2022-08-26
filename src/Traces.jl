# Traces wrapper for julia.

"Julia wrapper for the Traces C library."
module Traces


using Libdl: dlext

function depsdir(pkg::AbstractString)
    pkgdir = Base.find_package(pkg)
    pkgdir = abspath(joinpath(dirname(pkgdir), "..", "deps"))
    return pkgdir
end

const LIB_FILE = joinpath(depsdir("Traces"), "mintraceswrap." * dlext) 

const Nboolean = Cint

mutable struct TracesOptions
    getcanon::Nboolean        # make canong and canonlab?
    writeautoms::Nboolean     # write automorphisms?
    cartesian::Nboolean       # use cartesian rep for writing automs?
    digraph::Nboolean         # multiple edges or loops?
    defaultptn::Nboolean      # set lab,ptn,active for single cell?
    linelength::Cint          # max chars/line (excl. '\n') for output.
    outfile::Ptr{Cvoid}       # file for output, if any. FILE *outfile,
    strategy::Cint            # Only the value 0 is supported in this version.
    verbosity::Cint           # A level of verbosity of messages while Traces is running.
    generators::Ptr{Cvoid}    # This can be used to provide known automorphisms to Traces.
    userautomproc::Ptr{Cvoid} # void (*userautomproc)                                     # procedure called for each automorphism
                              # (int,int*,int*,int,int,int);
    usercanonproc::Ptr{Cvoid} # Cint  (*usercanonproc)                                    # procedure called for better labellings
                              # (graph*,int*,graph*,int,int,int,int);
    weighted::Nboolean        #  wegihted graph ? not used in this version 
end


const TracesOptions() = ccall((:defaultoptions_traces, LIB_FILE), TracesOptions, ())

const DEFAULTOPTIONS_TRACES = TracesOptions()

mutable struct TracesStats
    grpsize1::Cdouble        # /* size of group is */
    grpsize2::Cint           # /* grpsize1 * 10^grpsize2 */
    numgenerators::Cint      # /* number of generators found */
    numorbits::Cint          # /* number of orbits in group */
    treedepth::Cint          # /* The depth of the search tree */
    canupdates::Culong       # /* number of updates of best label */
    errstatus::Cint          # /* if non-zero : an error code */
    numnodes::Culong         # /* total number of nodes */
    interrupted::Culong      # /* The number of refinement operations aborted early. */
    peaknodes::Culong        # /* The maximum number of tree nodes simultaneously existing */

    function TracesStats()
        stats = new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        return stats
    end
end


mutable struct SparseGraph
    nde::Csize_t             # /* Number of directed edges (loops contribute only 1) */
    v::Ptr{Csize_t}          # /* Array of indexes into e[*] */
    nv::Cint                 # /* Number of vertices */
    d::Ptr{Cint}             # /* Array with out-degree of each vertex */
    e::Ptr{Cint}             # /* Array to hold lists of neighbours */
    w::Ptr{Cvoid}            # /* Not implemented, should be NULL. */
    vlen::Csize_t            # /* Sizes of arrays in units of type */
    dlen::Csize_t 
    elen::Csize_t 
    wlen::Csize_t
    
    function SparseGraph()
        sparsegraph = new(0, C_NULL, 0, C_NULL, C_NULL, C_NULL, 0, 0, 0, 0)
        return sparsegraph
    end
end





# Interface:

"""
    canong::SparseGraph

The canonical graph of the class of isomorphs that g is in. Only meaningful if
options.getcanon was 1

    labelling::Array{Cint}

if options.getcanon = 1, then this is the vertices of g in the order that they
should be relabelled to give canonical_graph.

    partition::Array{Cint}

colouring information for labels

    orbits::Array{Cint}

Orbits of the automorphism group

    stats::statsblk

stats related to the Traces run
"""
struct tracesreturn
    canong::SparseGraph
    labels::Array{Cint}
    partition::Array{Cint}
    orbits::Array{Cint}
    stats::TracesStats
end



function traces(g::SparseGraph,
                options = DEFAULTOPTIONS_TRACES::TracesOptions,
                labelling = nothing::Union{Cvoid, Array{Cint}},
                partition = nothing::Union{Cvoid, Array{Cint}})

    stats = TracesStats()

    # labelling and partition must be defined if defaultptn is not set and must not be defined if they are.
    @assert (labelling == nothing) == (options.defaultptn == 1)
    @assert (partition == nothing) == (options.defaultptn == 1)

    # Create some empty arrays for traces
    if options.defaultptn == 1
        labelling = zeros(Cint, g.nv)
        partition = zero(labelling)
    end

    # These don't need to be zero'd, I'm just doing it for debugging reasons.
    outgraph = SparseGraph()
    orbits = zero(labelling)

 
    ccall((:Traces, LIB_FILE), Cvoid,
          (Ref{SparseGraph}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{TracesOptions}, Ref{TracesStats}, Ref{SparseGraph}), g, labelling, partition, orbits, options, stats, outgraph)

    # Return everything nauty gives us.
    return tracesreturn(outgraph, labelling, partition, orbits, stats)
end
            
# helper methods:

using MetaGraphs, SimpleGraphs, Graphs

function fadjlist(g::GraphType) where GraphType <: Union{SimpleGraphs.AbstractSimpleGraph, Graphs.AbstractSimpleGraph}
    return g.fadjlist
end


function fadjlist(g::GraphType) where GraphType <: MetaGraphs.AbstractMetaGraph
    return g.graph.fadjlist
end

function to_sparse(g::GraphType) where GraphType <: Union{SimpleGraphs.AbstractSimpleGraph, Graphs.AbstractSimpleGraph, MetaGraphs.AbstractMetaGraph}
    num_vertices = nv(g)
    num_edges = ne(g)

    sparsegraph = SparseGraph()

    d = Array{Cint}(undef, num_vertices)
    v = Array{Csize_t}(undef, num_vertices)

    vi = 0
    vlen = 0
    dlen = 0
    for (i, neighbors) in enumerate(fadjlist(g))
        degree = length(neighbors)
        d[i] = convert(Cint, degree)
        dlen += 1
        v[i] = convert(Csize_t, vi)
        vlen += 1
        vi += degree
    end

    nde = 0
    for i in  1:num_vertices 
        nde =  nde > (v[i] + d[i]) ?  nde : (v[i] + d[i]) 
    end
    e = Array{Cint}(undef, nde)

 
    elen = 0 
    for (i, neighbors) in enumerate(fadjlist(g))
        vi = v[i] + 1
        for j in neighbors
            e[vi] = j - 1
            vi += 1
            elen += 1
        end
    end

    sparsegraph.nde =  convert(Csize_t, nde)
    sparsegraph.nv = convert(Csize_t, num_vertices)
    sparsegraph.d = pointer(d)
    sparsegraph.v = pointer(v)
    sparsegraph.e = pointer(e)
    sparsegraph.w = C_NULL
    sparsegraph.wlen = convert(Csize_t, 0)
    sparsegraph.vlen = convert(Csize_t, vlen)
    sparsegraph.dlen = convert(Csize_t, dlen)
    sparsegraph.elen = convert(Csize_t, elen)
    

    return sparsegraph
end

# return a dictionary of sets (orbits) of a SparseGraph 
function orbits(sparsegraph::SparseGraph, 
                labelling = nothing::Union{Array{Cint}, Nothing},
                partition = nothing::Union{Array{Cint}, Nothing})

    options =  TracesOptions()
    if labelling != nothing && partition != nothing
        options.defaultptn = 0
    end
    tracesreturn = traces(sparsegraph, options, labelling, partition)
    orbits = tracesreturn.orbits
    num_vertices = sparsegraph.nv

    orbit_class = Dict{Int, Set{Int}}()
    for v_ind in 1:num_vertices
        class_ind = orbits[v_ind] + 1 # C to julia index
        if !haskey(orbit_class, class_ind)
            orbit_class[class_ind] = Set{Int}()
            @assert(class_ind == v_ind)
        end
        push!(orbit_class[class_ind], v_ind) 
    end

    return orbit_class
end


# return array-like labelling, parition of a dictionary of sets (labels)
function to_label(labels::Dict{Int, Set{Int}})
    num_vertices = 0
    for label in values(labels)
        for v_ind in label
            num_vertices = num_vertices > v_ind  ? num_vertices : v_ind
        end
    end
    labelling = Array{Cint}(undef, num_vertices)
    partition = Array{Cint}(undef, num_vertices)
    ind = 1
    for label in values(labels)
        if length(label) != 0
            for v_ind in label
                labelling[ind] = convert(Cint, v_ind - 1) # julia to C index
                partition[ind] = convert(Cint, 1) 
                ind += 1
            end
            partition[ind - 1] = 0 
        end
    end
    return labelling, partition
end

end
