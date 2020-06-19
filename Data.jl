using SimpleWeightedGraphs

struct Data
    n::Int
    k::Int
    G::SimpleWeightedGraph
    degree::Array{Int}
    M::Int
    C::Float64
    instance::String
end

function Data(n, k, G, instance)
    degree = [ sum(G.weights[i, :]) for i = 1:n ]
    M = .5*sum(degree)
    C = M*(log(2.0*M) - 1)
    return Data(n, k, G, degree, M, C, instance)
end