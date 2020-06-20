using SimpleWeightedGraphs

struct Data
    # Number of samples
    n::Int
    
    # Number of communities
    k::Int

    # Graph
    G::SimpleWeightedGraph
    
    # Samples degree
    degree::Array{Int}
    
    # Number of edges
    M::Int

    # SBM constant
    C::Float64

    # Dataset name
    instance::String
end

function Data(n, k, G, instance)
    degree = [ sum(G.weights[i, :]) for i = 1:n ]
    M = .5*sum(degree)
    C = M*(log(2.0*M) - 1)
    return Data(n, k, G, degree, M, C, instance)
end