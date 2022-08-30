import Graphs
import Traces


function test_raw()
    sparse_g = Traces.SparseGraph()
    sparse_g.nv = 3
    sparse_g.nde = 4
    d = Array{Cint}([2,1,1])
    e = Array{Cint}([1,2,0,0,0,0,0])
    v = Array{Csize_t}([0, 3, 5])
    sparse_g.d = pointer(d)
    sparse_g.e = pointer(e)
    sparse_g.v = pointer(v)
    sparse_g.dlen = convert(Csize_t, 5)
    sparse_g.vlen = convert(Csize_t, 5)
    sparse_g.elen = convert(Csize_t, 7)

    print(sparse_g, "\n")
    print(d, "\n", e, "\n", v, "\n")
    
    orbit_class = Traces.orbits(sparse_g)

    print(orbit_class, "\n")    
end

function test_simpleline(n)
    g = Graphs.Graph()
    Graphs.add_vertex!(g)
    for i in  2:n
        Graphs.add_vertex!(g)
        Graphs.add_edge!(g, i-1, i)
    end

    sparse_g = Traces.to_sparse(g)

    orbit_class = Traces.orbits(sparse_g)

    print(orbit_class, "\n")
end

function test_cluster(n, c)
    g = Graphs.Graph()
    for i in  1:n
        Graphs.add_vertex!(g)
    end

    for i in 1:c
        for j in (i+1):c 
            Graphs.add_edge!(g, i, j)
        end
    end

    sparse_g = Traces.to_sparse(g)

    orbit_class = Traces.orbits(sparse_g)

    print(orbit_class, "\n")
end

function test_simple_graph()
    test_simpleline(10)
    test_simpleline(64)
    test_simpleline(122)
    test_cluster(10, 5)
    test_cluster(100, 30)
end

test_raw()
test_simple_graph()
