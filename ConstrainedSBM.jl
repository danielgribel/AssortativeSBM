using Convex
using SCS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
using Discreet

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

function compute_m(data, y)
    G = data.G
    k = data.k
    m = zeros(Int, k, k)
    for e in collect(edges(G))
        a = e.src
        b = e.dst
        m[ y[a], y[b] ] += e.weight
        m[ y[b], y[a] ] += e.weight
    end
    return m
end

function compute_kappa(data, y)
    G = data.G
    k = data.k
    kappa = zeros(Int, k)
    for i = 1:data.n
        kappa[ y[i] ] += sum(G.weights[i,:])
    end
    return kappa
end

function Solution(data, y)
    z = zeros(Int, data.n, data.k)
    # m = zeros(Int, data.k, data.k)
    # kappa = zeros(Int, data.k)
    
    m = compute_m(data, y)
    kappa = compute_kappa(data, y)

    for i = 1:data.n
        z[i, y[i]] = 1
    end

    # for e in collect(edges(G))
    #     a = e.src
    #     b = e.dst
    #     m[ y[a], y[b] ] += e.weight
    #     m[ y[b], y[a] ] += e.weight
    # end

    # for i = 1:data.n
    #     kappa[ y[i] ] += sum(data.G.weights[i,:])
    # end

    ll = 0.
    if is_assortative(m, kappa)
        println("***************** Assortative = YES")
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
        ll = 0.5*ll - C
    else
        println("***************** Assortative = NO")
        ll = -solveSBM(y, z, data)
    end
    return Solution(ll, data, y, z, m, kappa)
end

function solveSBM(y, z, data)
    n = data.n
    k = data.k
    A = data.A

    w = Variable(k, k) # SBM parameters
    
    # w.value = rand(k, k)

    # println("-- W declared.")

    # f = 0.5*sum(
    #     sum( ( A[i,j]*log(w[r,s]) - Q[i,j]*w[r,s] )*z[i,r]*z[j,s] for i=1:n, j=1:n )
    #     for r=1:k, s=1:k 
    # )

    f  = 0.5*sum( ( A[i,i]*log(w[ y[i], y[i] ]) - Q[i,i]*w[ y[i], y[i] ] ) for i=1:n )
    f += sum( ( A[i,j]*log(w[ y[i], y[j] ]) - Q[i,j]*w[ y[i], y[j] ] ) for i=1:n for j=(i+1):n )

    # println("-- Model built.")

    problem = maximize(f)

    # println("-- Problem declared.")

    # a = Array((1:k))
    # for r in a
    #     for s in a[a .!= r]
    #         problem.constraints += [w[r,r] > w[r,s]]
    #     end
    # end
    
    problem.constraints += w[1,1] >= (w[1,2] + Tol)
    problem.constraints += w[1,1] >= (w[2,1] + Tol)
    problem.constraints += w[2,2] >= (w[1,2] + Tol)
    problem.constraints += w[2,2] >= (w[2,1] + Tol)
    
    problem.constraints += w[1,1] >= Tol
    problem.constraints += w[1,2] >= Tol
    problem.constraints += w[2,1] >= Tol
    problem.constraints += w[2,2] >= Tol

    problem.constraints += w[1,2] == w[2,1]

    # problem.constraints += w[1,2] <= 0.5
    # problem.constraints += w[2,1] <= 0.5

    # println("-- Constraints added.")

    solve!(problem, SCSSolver(verbose = false))

    # solve!(problem, SCSSolver(
    #     verbose = false,
    #     # normalize = 0,    # boolean, heuristic data rescaling: 1
    #     # scale = 1.0,      # if normalized, rescales by this factor: 5
    #     #rho_x = 1e-3,      # x equality constraint scaling: 1e-3
    #     max_iters = 5,      # maximum iterations to take: 2500
    #     eps = 1e-4,         # convergence tolerance: 1e-3
    #     # alpha = 1.8,      # relaxation parameter: 1.8
    #     # cg_rate = 1
    # ))
    
    # println("-- Problem solved.")
    println("Optval = ", problem.optval)
    println("Omega = ", w.value)
    return problem.optval
end

function evalrelocate(sol, i, t)
    # source cluster
    src = sol.y[i]

    m_ = copy(sol.m)

    for v in neighbors(sol.data.G, i)
        m_[sol.y[v], src] -= sol.data.G.weights[i, v]
        m_[src, sol.y[v]] -= sol.data.G.weights[i, v]
        if i == v
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
    sbm_cost = 0.5*sbm_cost - C
        
    if sbm_cost < (sol.ll + Tol)
        if is_assortative(m_, kappa_)
            # update likelihood
            sol.ll = sbm_cost
            # update parameters
            sol.m = copy(m_)
            # update assignments
            sol.z[i, src] = 0
            sol.z[i, t] = 1
            sol.y[i] = t
            sol.kappa = copy(kappa_)
            # println("--relocate ", i, " ", t)
        else
            y_ = copy(sol.y)
            z_ = copy(sol.z)
            y_[i] = t
            z_[i, src] = 0
            z_[i, t] = 1
            sbm_cost = -solveSBM(y_, z_, sol.data)
            if sbm_cost < (sol.ll + Tol)
                # update likelihood
                sol.ll = sbm_cost
                # update parameters
                sol.m = copy(m_)
                # update assignments
                sol.z[i, src] = 0
                sol.z[i, t] = 1
                sol.y[i] = t
                sol.kappa = copy(kappa_)
                # println("--relocate ", i, " ", t)
            end
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
    counter = 0

    while ((prev_ll - curr_ll) > Tol) && it_changed
        rdm_samples = randperm(tw, n)
        it_changed = false
        for i in rdm_samples
            prev = y[i]
            rdm_clusters = randperm(tw, k)
            for c in rdm_clusters
                if y[i] != c && prev != c
                    evalrelocate(sol, i, c)
                end
                if sol.ll < (curr_ll + Tol)
                    prev_ll = curr_ll
                    curr_ll = sol.ll
                    it_changed = true
                end
            end
        end
        counter += 1
    end
    println(counter)
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
    ll = 0.5*ll - C
    return ll
end

function get_omega(m, kappa)
    k = length(kappa)
    w = [ 2.0*M*(m[r, s])/(kappa[r] * kappa[s]) for r=1:k, s=1:k ]
    return w
end

function is_assortative(m, kappa)
    return is_strongly_assortative(m, kappa)
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

function run()
    n = data.n
    k = data.k
    A = data.A
    max_it = 10
    lls = zeros(Float64, max_it)

    for i=1:max_it
        # create initial solution
        y = rand([1,k], n)
        m = compute_m(data, y)
        kappa = compute_kappa(data, y)
        while !is_assortative(m, kappa)
            y = rand([1,k], n)
            m = compute_m(data, y)
            kappa = compute_kappa(data, y)
        end
        sol = Solution(data, y)
        ll_initial = sol.ll
        localsearch(sol)
        println("Likelihood (", i, "): ", ll_initial, " >> ", sol.ll)
        nmi = mutual_information(sol.y, label; normalize = true)
        println("NMI: ", nmi)
        w = get_omega(sol.m, sol.kappa)
        println("W: ", w)
        lls[i] = sol.ll
    end
    println("minimum: ", minimum(lls))
end


### MAIN

INPUT_FILE = "data/Data-2-8-100-1000.link"

LABEL_FILE = "data/Data-2-8-100-1000.label"

# tolerance epsilon
Tol = 1e-4

# number of samples
n = 100

# data dimensionality
d = 8

# number of clusters
k = 2

# graph (adjacency matrix)
A = zeros(Int, n, n)

nb_pairs = countlines(INPUT_FILE)

G = SimpleWeightedGraph(n)

label = zeros(Int, n)

open(INPUT_FILE) do file
    for ln in eachline(file)
        lnsplit = split(ln, " ")
        a = parse(Int, lnsplit[1])
        b = parse(Int, lnsplit[2])
        e = parse(Int, lnsplit[3])
        A[a, b] = e
        A[b, a] = e
        add_edge!(G, a, b, e)
    end
end

open(LABEL_FILE) do file
    ct = 1
    for ln in eachline(file)
        label[ct] = parse(Int, ln)
        ct += 1
    end
end

data = Data(n, d, k, A, G)

# total number of edges
M = .5*sum(A)

# constant
C = M*(log(2*M) - 1)

# (k_i k_j)/2M matrix
Q = ( sum(A[i,:] for i in 1:n) * transpose(sum(A[j,:] for j in 1:n)) )/(2*M)

lowerbound = 0.5*sum(Q)

run()

# y = [1,2,2,1,2,1,1,1,2,2,1,2,2,1,1,1,2,2,2,1,2,1,1,1,1,2,1,1,2,2,2,2,1,1,2,2,1,2,2,1,1,2,2,1,2,1,1,1,2,1,2,1,2,1,2,2,2,1,2,1,2,2,1,1,2,2,2,2,2,1,2,1,2,2,2,2,2,2,1,1,1,1,2,1,2,1,1,1,2,2,2,2,2,2,1,2,2,2,1,1,2,2,2,2,1,1,1,2,1,2,1,2,2,1,1,1,1,1,1,2,1,1,2,2,1,1,1,2,2,1,1,1,2,1,2,2,2,1,2,2,2,2,1,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,2,1,2,2,2,1,1,2,2,2,2,2,2,2,1,2,2,1,1,2,2,1,1,2,2,2,1,2,2,1,1,1,2,2,1,1,1]
