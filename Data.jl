using SimpleWeightedGraphs

struct Data
    n::Int
    d::Int
    k::Int
    G::SimpleWeightedGraph
    degree::Array{Int}
end

function Data(n, d, k, G)
    degree = [ sum(G.weights[i,:]) for i=1:n ]
    return Data(n, d, k, G, degree)
end