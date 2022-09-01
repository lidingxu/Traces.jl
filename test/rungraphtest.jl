import Graphs
import Traces




function test_simpleline(n)
    g = Graphs.Graph()
    Graphs.add_vertex!(g)
    for i in  2:n
        Graphs.add_vertex!(g)
        Graphs.add_edge!(g, i-1, i)
    end

    _, automorphism, orbit_class = Traces.traces(g)

    print("\n orbit:",orbit_class, "\n")
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

    _, automorphism, orbit_class = Traces.traces(g, false, true)

    print("\n orbit:", orbit_class, "\n autom:", automorphism , "\n")
end

function test_simple_graph()
    test_simpleline(10)
    test_simpleline(64)
    test_simpleline(122)
    test_cluster(10, 5)
    test_cluster(100, 30)
end

test_simple_graph()
