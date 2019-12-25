using Convex
using SCS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs

struct Data
    n::Int
    d::Int
    k::Int
    A::Array{Int}
    G::SimpleWeightedGraph
end

mutable struct Solution
    ll::Float64
    data::Data
    y::Array{Int}
    z::Array{Int}
    m::Array{Int}
    kappa::Array{Int}
end

function Data(n, d, k, A, G)
    return Data(n, d, k, A, G)
end

function Solution(data, y)
    z = zeros(Int, data.n, data.k)
    m = zeros(Int, data.k, data.k)
    kappa = zeros(Int, data.k)
    
    for i = 1:data.n
        z[i, y[i]] = 1
    end

    for e in collect(edges(G))
        a = e.src
        b = e.dst
        m[ y[a], y[b] ] += e.weight
        m[ y[b], y[a] ] += e.weight
    end

    for r = 1:data.k
        for s = 1:data.k
            kappa[r] += m[r, s]
        end
    end

    ll = 0.

    for r = 1:data.k
        for s = r:data.k
            if (m[r,s] == 0) || (kappa[r]*kappa[s] == 0)
               contribM = 0.
            else
                denM = kappa[r] * kappa[s]
                contribM = m[r,s] * log(float(m[r,s])/denM)
            end
            ll -= contribM
            if r != s
                ll -= contribM
            end
        end
    end
    ll *= 0.5
    ll -= C

    return Solution(ll, data, y, z, m, kappa)
end

function solveSBM(y, data)
    n = data.n
    k = data.k

    # number of edges
    m = .5 * sum(data.A)
    
    kmatrix = ( sum(data.A[i,:] for i in 1:n) * transpose(sum(data.A[j,:] for j in 1:n)) )/(2*m)

    z = Array{Int, 2}(undef, n, k)
    fill!(z, 0)

    for i in 1:n
        z[i, y[i]] = 1
    end

    w = Variable(k, k)

    obj = .5*sum(
        sum(( data.A[i,j]*log(w[r,s]) - w[r,s]*kmatrix[i,j] )*z[i,r]*z[j,s] for i in 1:n, j in 1:n )
        for r in 1:k, s in 1:k 
    )

    problem = maximize(obj)

    a = Array((1:k))

    for r in a
        for s in a[a .!= r]
            problem.constraints += [w[r,r] > w[r,s]]
        end
    end

    solve!(problem, SCSSolver())
    return problem.optval
end

function evalrelocate(sol, i, t)
    # source center
    src = sol.y[i]

    m_ = copy(sol.m)

    for v in neighbors(sol.data.G, i)
        m_[sol.y[v], src] -= sol.data.G.weights[i, v]
        m_[src, sol.y[v]] -= sol.data.G.weights[i, v]
        if(i == v)
            m_[t, t] += 2.0*sol.data.G.weights[i, i]
        else
            m_[t, sol.y[v]] += sol.data.G.weights[i, v]
            m_[sol.y[v], t] += sol.data.G.weights[i, v]
        end
    end

    sbm_cost = 0.

    kappa_ = copy(sol.kappa)
    kappa_[ src ] -= sum(sol.data.G.weights[i,:]) 
    kappa_[ t ] += sum(sol.data.G.weights[i,:])

    # SBM likelihood
    for r = 1:sol.data.k
        for s = r:sol.data.k
            if (m_[r,s] == 0) || (kappa_[r]*kappa_[s] == 0)
               contribM = 0.
            else
                denM = kappa_[r] * kappa_[s]
                contribM = m_[r,s] * log(float(m_[r,s])/denM)
            end
            sbm_cost -= contribM
            if r != s
                sbm_cost -= contribM
            end
        end
    end
    
    sbm_cost *= 0.5
    sbm_cost -= C

    if sbm_cost < (sol.ll + Tol) && is_strongly_assortative(m_, kappa_)
        # update likelihood
        sol.ll = sbm_cost

        # update parameters
        sol.m = copy(m_)

        # update assignments
        sol.z[i, src] = 0
        sol.z[i, t] = 1
        sol.y[i] = t
        sol.kappa = copy(kappa_)
    end

    if sbm_cost < (sol.ll + Tol) && !is_strongly_assortative(m_, kappa_)
        y_ = copy(sol.y)
        y_[i] = t
        cost_constrained = -solveSBM(y_, sol.data)
        if cost_constrained < (sol.ll + Tol)
            # update likelihood
            sol.ll = cost_constrained

            # update parameters
            sol.m = copy(m_)

            # update assignments
            sol.z[i, src] = 0
            sol.z[i, t] = 1
            sol.y[i] = t
            sol.kappa = copy(kappa_)
        end
    end
end

function localsearch(sol)
    y = sol.y
    n = sol.data.n
    k = sol.data.k
    tw = MersenneTwister(1234)
    prev_ll = Inf
    curr_ll = sol.ll
    it_changed = true

    while (prev_ll - curr_ll) > Tol && it_changed
        rdm_samples = randperm(tw, n)
        it_changed = false
        for i in rdm_samples
            prev_c = y[i]
            rdm_clusters = randperm(tw, k)
            for c in rdm_clusters
                if y[i] != c && prev_c != c
                    evalrelocate(sol, i, c)
                end
                if sol.ll < (curr_ll + Tol)
                    prev_ll = curr_ll
                    curr_ll = sol.ll
                    it_changed = true
                end
            end
        end
    end
end

function likelihood(sol)
    ll = 0.
    for r = 1:sol.data.k
        for s = r:sol.data.k
            if (sol.m[r,s] == 0) || (sol.kappa[r]*sol.kappa[s] == 0)
               contribM = 0.
            else
                denM = sol.kappa[r] * sol.kappa[s]
                contribM = sol.m[r,s] * log(float(sol.m[r,s])/denM)
            end
            ll -= contribM
            if r != s
                ll -= contribM
            end
        end
    end
    ll *= 0.5
    ll -= C
    return ll
end

function get_omega(m, kappa)
    k = length(kappa)
    w = [ (m[r, s])/(kappa[r] * kappa[s]) for r in 1:k, s in 1:k ]
    return w
end

function is_strongly_assortative(m, kappa)
    w = get_omega(m, kappa)    
    if minimum(diag(w)) > maximum(w - Diagonal(w))
        return true
    end
    return false
end

function is_weakly_assortative(m, kappa)
    k = length(kappa)
    w = get_omega(m, kappa)
    a = Array((1:k))
    for r = 1:k
        if w[r, r] <= w[r, a[a .!= r]]
            return false
        end
    end
    return true
end



INPUT_FILE = "data/Sample.link"

# tolerance epsilon
Tol = 1e-4

# number of samples
n = 8

# data dimensionality
d = 2

# number of clusters
k = 2

# graph (adjacency matrix)
A = zeros(Int, n, n)

nb_pairs = countlines(INPUT_FILE)

G = SimpleWeightedGraph(nb_pairs)

open(INPUT_FILE) do file
    for ln in eachline(file)
        lnsplit = split(ln, " ")
        a = parse(Int8, lnsplit[1])
        b = parse(Int8, lnsplit[2])
        e = parse(Int8, lnsplit[3])
        A[a, b] = e
        A[b, a] = e
        add_edge!(G, a, b, e)
    end
end

data = Data(n, d, k, A, G)

# total number of edges
m = .5*sum(A)

# constant
C = m*(log(2*m) - 1)

# (k_i k_j)/2m matrix
Q = ( sum(A[i,:] for i in 1:n) * transpose(sum(A[j,:] for j in 1:n)) )/(2*m)

y = zeros(Int, data.n)

y = [1, 2, 1, 1, 1, 2, 2, 2]

sol = Solution(data, y)

println(sol.ll)

localsearch(sol)

println(sol.ll)

println(sol.y)