using Convex
using SCS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
using Discreet
using StatsBase

struct Data
    n::Int
    d::Int
    k::Int
    A::Array{Int}
    G::SimpleWeightedGraph
    degree::Array{Int}
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
    degree = [ sum(A[i,:]) for i=1:n ]
    return Data(n, d, k, A, G, degree)
end

function compute_m(data, y)
    G = data.G
    A = data.A
    k = data.k
    m = zeros(Int, k, k)
    for i=1:data.n
        for j=1:data.n
            m[ y[i], y[j] ] += A[i,j]
        end
    end
    return m
end

function compute_kappa(data, y)
    G = data.G
    k = data.k
    kappa = zeros(Int, k)
    for i = 1:data.n
        kappa[ y[i] ] += data.degree[i]
    end
    return kappa
end

function Solution(data, y)
    z = zeros(Int, data.n, data.k)
    m = compute_m(data, y)
    kappa = compute_kappa(data, y)

    for i = 1:data.n
        z[i, y[i]] = 1
    end

    ll = 0.
    if is_assortative(m, kappa)
        println("***************** Assortative = YES")
        for r = 1:data.k
            for s = r:data.k
                if (m[r,s] == 0) || (kappa[r]*kappa[s] == 0)
                contribM = 0.
                else
                    denM = kappa[r] * kappa[s]
                    contribM = m[r,s] * log(m[r,s]/denM)
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
        ll = -solveSBM(m, kappa, data)
    end
    return Solution(ll, data, y, z, m, kappa)
end

function solveSBM(m_, kappa_, data)
    n = data.n
    k = data.k
    A = data.A
    m = m_
    kap = zeros(Float64, k, k)
    
    for r=1:k
        for s=r:k
            if m_[r,s] == 0 || kappa_[r] == 0
                return -Inf
            end
        end
    end

    for r=1:k
        for s=r:k
            kap[r,s] = (kappa_[r]*kappa_[s])/(2*M)
            kap[s,r] = kap[r,s]
        end
    end

    # SBM parameters
    w = Variable(k, k)

    # Objective
    f  = 0.5*sum( ( m[r,r]*log(w[r,r]) - w[r,r]*kap[r,r] ) for r=1:k )
    f += sum( ( m[r,s]*log(w[r,s]) - w[r,s]*kap[r,s] ) for r=1:k for s=(r+1):k )

    problem = maximize(f)

    # Model constraints
    for r = 1:k
        for s = r:k
            problem.constraints += [ w[r,s] >= Tol ]
            if r != s
                problem.constraints += [ w[r,s] == w[s,r] ]
            end
        end
    end

    for q = 1:k
        for r = 1:k
            for s = (r+1):k
                problem.constraints += [ w[q,q] >= (w[r,s] + Tol) ]
            end
        end
    end

    solve!(problem, SCSSolver(verbose = false))
    # println("Optval = ", problem.optval)
    # println("Omega = ", w.value)
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
    kappa_[ src ] -= sol.data.degree[i]
    kappa_[ t ] += sol.data.degree[i]

    for r = 1:sol.data.k
        for s = r:sol.data.k
            if (m_[r,s] == 0) || (kappa_[r]*kappa_[s] == 0)
            contribM = 0.
            else
                denM = kappa_[r] * kappa_[s]
                contribM = m_[r,s] * log(m_[r,s]/denM)
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
            sbm_cost = -solveSBM(m_, kappa_, sol.data)
            if sbm_cost < (sol.ll + Tol)
                # update likelihood
                sol.ll = sbm_cost
                # update parameters
                sol.m = copy(m_)
                # update assignments
                sol.y[i] = t
                sol.z[i, src] = 0
                sol.z[i, t] = 1
                sol.kappa = copy(kappa_)
                # println("------------Relocate after SOLVE_SBM ", i, " ", t)
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
    # println(counter)
end

function likelihood(sol)
    ll = 0.
    for r = 1:sol.data.k
        for s = r:sol.data.k
            if (sol.m[r,s] == 0) || (sol.kappa[r]*sol.kappa[s] == 0)
               contribM = 0.
            else
                denM = sol.kappa[r] * sol.kappa[s]
                contribM = sol.m[r,s] * log(sol.m[r,s]/denM)
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
    for r = 1:k
        for s = 1:k
            if r != s && w[r,r] <= w[r,s]
                return false
            end
        end
    end
    return true
end

function run()
    n = data.n
    k = data.k
    A = data.A
    max_it = 10
    best_nmi = 0.0
    best_ll = Inf

    for i=1:max_it
        # create initial solution
        y = sample(1:k, n)
        m = compute_m(data, y)
        kappa = compute_kappa(data, y)
        # while !is_assortative(m, kappa)
        #     y = sample(1:k, n)
        #     m = compute_m(data, y)
        #     kappa = compute_kappa(data, y)
        # end
        sol = Solution(data, y)
        ll_initial = sol.ll
        println(ll_initial)
        localsearch(sol)
        nmi = mutual_information(sol.y, label; normalize = true)
        w = get_omega(sol.m, sol.kappa)
        println("Likelihood (", i, "): ", ll_initial, " >> ", sol.ll)
        println("NMI: ", nmi)
        println("W: ", w)
        if sol.ll < best_ll
            best_ll = sol.ll
            best_nmi = nmi
        end
    end
    println("minimum: ", best_ll, " ", best_nmi)
end


### MAIN

INPUT_FILE = "data/Data-4-8-200-1001.link"

LABEL_FILE = "data/Data-4-8-200-1001.label"

# tolerance epsilon
Tol = 1e-4

# number of samples
n = 200

# data dimensionality
d = 8

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

# number of clusters
k = length(unique(collect(label)))

data = Data(n, d, k, A, G)

# total number of edges
M = .5*sum(A)

# SBM constant
C = M*(log(2*M) - 1)

# (k_i k_j)/2M constant matrix
Q = (data.degree * transpose(data.degree))/(2*M)

run()