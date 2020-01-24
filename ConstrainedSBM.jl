using Convex
using SCS
using ECOS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
using Discreet
using StatsBase

include("Data.jl")

mutable struct Solution
    ll::Float64
    data::Data
    y::Array{Int}
    z::Array{Int}
    m::Array{Int}
    kappa::Array{Int}
end

function Solution(data, y)
    z = zeros(Int, data.n, data.k)
    m = compute_m(data, y)
    kappa = compute_kappa(data, y)
    
    for i = 1:data.n
        z[i, y[i]] = 1
    end

    # q = [ sum(z[:,r]) for r = 1:data.k ]

    if is_assortative(m, kappa)
        println("***************** Assortative = YES")
        ll = ll_DCSBM(data, m, kappa)
    else
        println("***************** Assortative = NO")
        ll = solve_DCSBM(data, m, kappa)
    end
    return Solution(ll, data, y, z, m, kappa)
end

function compute_m(data, y)
    A = data.A
    m = zeros(Int, data.k, data.k)
    for i=1:data.n
        m[ y[i], y[i] ] += A[i,i]
        for j=(i+1):data.n
            m[ y[i], y[j] ] += A[i,j]
            m[ y[j], y[i] ] += A[j,i]
        end
    end
    return m
end

function compute_kappa(data, y)
    kappa = zeros(Int, data.k)
    [ kappa[ y[i] ] += data.degree[i] for i=1:data.n ]
    return kappa
end

function ll_SBM(data, m, q)
    ll = 0.
    for r = 1:data.k
        for s = r:data.k
            ratio = m[r,s]/(q[r] * q[s])
            if ratio < 1e-6
                ratio = 1e-10
            end
            contrib = m[r,s] * log(ratio)
            ll -= contrib
            if r != s
                ll -= contrib
            end
        end
    end
    ll = 0.5*ll - C
    return ll
end


function ll_DCSBM(data, m, kappa)
    ll = 0.
    for r = 1:data.k
        for s = r:data.k
            ratio = m[r,s]/(kappa[r] * kappa[s])
            if ratio < 1e-6
                ratio = 1e-10
            end
            contrib = m[r,s] * log(ratio)
            ll -= contrib
            if r != s
                ll -= contrib
            end
        end
    end
    ll = 0.5*ll - C
    return ll
end

function solve_SBM(data, m, q)
    n = data.n
    k = data.k
    A = data.A
    
    # SBM parameters
    w = Variable(k, k)

    # Objective
    f  = 0.5*sum( ( m[r,r]*log(w[r,r]) - q[r]*q[r]*w[r,r] ) for r=1:k )
    f += sum( ( m[r,s]*log(w[r,s]) - q[r]*q[s]*w[r,s] ) for r=1:k for s=(r+1):k )

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

    for a = 1:k
        for r = 1:k
            for s = (r+1):k
                problem.constraints += [ w[a,a] >= (w[r,s] + Tol) ]
            end
        end
    end

    solve!(problem, ECOSSolver(verbose = false))
    # println("Optval = ", -problem.optval)
    # println("Omega = ", w.value)
    return -problem.optval
end

function add_feasibility_constraints(problem, w, k)
    for r = 1:k
        for s = r:k
            # w > 0 constraint
            problem.constraints += [ w[r,s] >= Tol ]
            # w[r,s] = w[s,r] constraint
            if r != s
                problem.constraints += [ w[r,s] == w[s,r] ]
            end
        end
    end
end

function add_strong_constraints(problem, w, k)
    for a = 1:k
        for r = 1:k
            for s = (r+1):k
                # w[a,a] > w[r,s] for all a, for all r,s; r != s
                problem.constraints += [ w[a,a] >= (w[r,s] + Tol) ]
            end
        end
    end
end

function add_weak_constraints(problem, w, k)
    for a = 1:k
        for r = 1:k
            if r != a
                # w[a,a] > w[a,r] for all a, for all r, r != a
                problem.constraints += [ w[a,a] >= (w[a,r] + Tol) ]
            end
        end
    end
end

function solve_DCSBM(data, m, kappa_)
    k = data.k
    kap = zeros(Float64, k, k)
    
    for r = 1:k
        for s = r:k
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

    # Add feasibility constraints: w[r,s] > 0 and w[r,s] == w[s,r]
    add_feasibility_constraints(problem, w, k)

    # Add model constraints
    add_strong_constraints(problem, w, k)

    # Solve the constrained convex problem
    solve!(problem, ECOSSolver(verbose = false))

    # println("Optval = ", -problem.optval)
    # println("Omega = ", w.value)
    
    return -problem.optval
end

function eval_relocate(sol, i, t)
    # source cluster
    src = sol.y[i]

    m_ = copy(sol.m)

    for v in neighbors(sol.data.G, i)
        e = Int(sol.data.G.weights[i, v])
        if i == v
            m_[src, src] -= e
            m_[t, t] += e
        else
            m_[sol.y[v], src] -= e
            m_[src, sol.y[v]] -= e
            m_[t, sol.y[v]] += e
            m_[sol.y[v], t] += e
        end
    end

    sbm_cost = 0.

    kappa_ = copy(sol.kappa)
    kappa_[ src ] -= sol.data.degree[i]
    kappa_[ t ] += sol.data.degree[i]

    sbm_cost = ll_DCSBM(sol.data, m_, kappa_)

    if sbm_cost < sol.ll
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
        else
            sbm_cost = solve_DCSBM(sol.data, m_, kappa_)
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
            end
        end
    end
end

function get_omega_SBM(m, q)
    k = length(q)
    w = [ (m[r, s])/(q[r] * q[s]) for r=1:k, s=1:k ]
    return w
end

function get_omega_DCSBM(m, kappa)
    k = length(kappa)
    w = [ 2.0*M*(m[r, s])/(kappa[r] * kappa[s]) for r=1:k, s=1:k ]
    return w
end

function is_assortative(m, kappa)
    if CONSTRAINED
        return is_strongly_assortative(m, kappa)
    end
    return true
end

function is_strongly_assortative(m, kappa)
    w = get_omega_DCSBM(m, kappa)
    if minimum(diag(w)) > maximum(w - Diagonal(w))
        return true
    end
    return false
end

function is_weakly_assortative(m, kappa)
    k = length(kappa)
    w = get_omega_DCSBM(m, kappa)
    for r = 1:k
        for s = 1:k
            if r != s && (w[r,r] <= w[r,s])
                return false
            end
        end
    end
    return true
end

function localsearch(sol)
    n = sol.data.n
    k = sol.data.k
    tw = MersenneTwister(1234)
    prev_ll = Inf
    curr_ll = sol.ll
    it_changed = true
    
    while ((prev_ll - curr_ll) > Tol) && it_changed
        rdm_samples = randperm(tw, n)
        it_changed = false
        for i in rdm_samples
            prev = sol.y[i]
            rdm_clusters = randperm(tw, k)
            for c in rdm_clusters
                if (sol.y[i] != c) && (prev != c)
                    eval_relocate(sol, i, c)
                end
                if sol.ll < (curr_ll - Tol)
                    prev_ll = curr_ll
                    curr_ll = sol.ll
                    it_changed = true
                end
            end
        end
    end
end

function run(max_it)
    best_nmi = 0.0
    best_ll = Inf
    best_solution = nothing

    ct_weakly = 0
    ct_strongly = 0

    for i=1:max_it
        # create initial solution
        y = sample(1:data.k, data.n)
        sol = Solution(data, y)
        ll_initial = sol.ll
        localsearch(sol)
        nmi = mutual_information(sol.y, label; normalize = true)
        println("Likelihood (", i, "): ", ll_initial, " >> ", sol.ll)
        println("NMI: ", nmi)
        
        w = get_omega_DCSBM(sol.m, sol.kappa)
        omega_ii = diag(w)
        omega_ij = w - Diagonal(w)
        [ println("w(", c, "): ", omega_ii[c], ", ", maximum(omega_ij[c,:])) for c = 1:k ]
        
        if sol.ll < best_ll
            best_ll = sol.ll
            best_nmi = nmi
            best_solution = sol
        end

        if is_weakly_assortative(sol.m, sol.kappa)
            ct_weakly = ct_weakly + 1
        end
        if is_strongly_assortative(sol.m, sol.kappa)
            ct_strongly = ct_strongly + 1
        end

    end

    best_w = get_omega_DCSBM(best_solution.m, best_solution.kappa)
    omega_ii = diag(best_w)
    omega_ij = best_w - Diagonal(best_w)
    yB = []
    for i=1:data.n
        if data.degree[i] != 0
            push!(yB, best_solution.y[i])
        end
    end
    nmiB = mutual_information(yB, labelB; normalize = true)

    println("----------- Best solution:")
    # [ println("w(", i, "): ", omega_ii[i], ",", maximum(omega_ij[i,:])) for i = 1:data.k ]
    [ println("w(", i, "): ", best_w[i,:]) for i = 1:data.k ]
    println("Size = ", [ sum(best_solution.z[:,r]) for r = 1:data.k ])
    println("Likelihood = ", best_ll)
    println("NMI = ", best_nmi, " ", nmiB)
    println("Assignment = ", best_solution.y)
    println("% assortative solutions = ", ct_weakly/max_it, " ", ct_strongly/max_it)
end


###### MAIN

CONSTRAINED = true

INSTANCE = "W-3-8-200-1001"

EDGES_FILE = "data/" * INSTANCE * ".link"
LABEL_FILE = "data/" * INSTANCE * ".label"

# tolerance epsilon
Tol = 1e-4

nb_pairs = countlines(EDGES_FILE)

label = Int[]

open(LABEL_FILE) do file
    for ln in eachline(file)
        push!(label, parse(Int, ln))
    end
end

# number of nodes
n = length(label)

# number of clusters
k = length(unique(collect(label)))

# graph (adjacency matrix)
A = zeros(Int, n, n)

G = SimpleWeightedGraph(n)

open(EDGES_FILE) do file
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

data = Data(n, 8, k, A, G)

# total number of edges
M = .5*sum(A)

# SBM constant
C = M*(log(2*M) - 1)

# (k_i k_j)/2M constant matrix for the DC-SBM
Q = (data.degree * transpose(data.degree))/(2*M)

labelB = []
for i=1:n
    if data.degree[i] != 0
        push!(labelB, label[i])
    end
end

run(500)

ground_truth = Solution(data, copy(label))

omega = get_omega_DCSBM(ground_truth.m, ground_truth.kappa)
omega_ii = diag(omega)
omega_ij = omega - Diagonal(omega)

for i=1:k
    println(omega_ii[i], ", ", maximum(omega_ij[i,:]))
end

println("ll = ", ll_DCSBM(data, ground_truth.m, ground_truth.kappa))