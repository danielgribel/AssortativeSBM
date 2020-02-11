using Convex
# using SCS
using ECOS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
using Discreet
using StatsBase
using Clustering
using CSV

import Base.iterate,Base.length

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
        ll = ll_DCSBM(data, m, kappa)
    else
        ll = solve_DCSBM(data, m, kappa)
    end
    return Solution(ll, data, y, z, m, kappa)
end

function compute_m(data, y)
    G = data.G
    m = zeros(Int, data.k, data.k)
    for i=1:data.n
        # m[ y[i], y[i] ] += A[i,i]
        m[ y[i], y[i] ] += Int(G.weights[i, i])
        for j=(i+1):data.n
            # m[ y[i], y[j] ] += A[i,j]
            # m[ y[j], y[i] ] += A[j,i]
            m[ y[i], y[j] ] += Int(G.weights[i, j])
            m[ y[j], y[i] ] += Int(G.weights[j, i])
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
        # return is_relaxed_assortative(m, kappa)
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

function is_relaxed_assortative(m, kappa)
    k = length(kappa)
    w = get_omega_DCSBM(m, kappa)
    count = 0
    for r = 1:k
        arr = findall(w[r, :] .== maximum(w[r, :]))
        if length(arr) == 1 && r == arr[1]
            count += 1
        end
    end
    if (count/k) >= PERC
        return true
    end
    return false
end

function localsearch(sol)
    n = sol.data.n
    k = sol.data.k
    prev_ll = Inf
    curr_ll = sol.ll
    it_changed = true
    
    while ((prev_ll - curr_ll) > Tol) && it_changed
        rdm_samples = randperm(Mt, n)
        it_changed = false
        for i in rdm_samples
            prev = sol.y[i]
            rdm_clusters = randperm(Mt, k)
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

function initial_assignment(data)
    rdm_order = randperm(Mt, data.n)
    y = zeros(Int, data.n)
    for i = 1:data.n
        y[ rdm_order[i] ] = ceil(data.k*i/data.n)
    end
    return y
end

function run(max_it)
    # best_nmi = 0.0
    # best_crand = -1.0
    # best_ll = Inf
    # best_solution = nothing

    for i=1:max_it
        # create initial solution
        y = initial_assignment(data)
        sol = Solution(data, y)
        ll_initial = sol.ll
        t1 = time_ns()
        localsearch(sol)
        elapsed_time = (time_ns() - t1)/1.0e9

        nmi = mutual_information(sol.y, label; normalize = true)
        crand = randindex(sol.y, label)[1]
        w = get_omega_DCSBM(sol.m, sol.kappa)
        omega_ii = diag(w)
        omega_ij = w - Diagonal(w)
        
        # if sol.ll < best_ll
        #     best_ll = sol.ll
        #     best_nmi = nmi
        #     best_crand = crand
        #     best_solution = sol
        # end

        count_csbm = 0
        for r = 1:data.k
            arr_csbm = findall(w[r, :] .== maximum(w[r, :]))
            arr_ground = findall(w_ground[r, :] .== maximum(w_ground[r, :]))
            if length(arr_csbm) == 1 && r == arr_csbm[1]
                count_csbm += 1
            end
        end

        line = INSTANCE * " "
        line *= string(SEED) * " "
        line *= string(data.k) * " "
        line *= string(sol.ll) * " "
        line *= string(nmi) * " "
        line *= string(count_csbm/data.k) * " "
        line *= string(elapsed_time) * " "
        for r = 1:data.k
            for s = 1:data.k
                line *= string(w[r,s]) * " "
            end
        end
        line *= "\n"
        write_output(line)
    end
end

function write_output(line)
    io = open(OUTPUT_FILE, "a")
    write(io, line)
    close(io)
end

struct Combinations{T}
    itr::Vector{T}
    count::Int64
    itrsize::Int64
    function Combinations(itr::Vector{T},count::Int) where T
        new{T}(itr,Int64(count),length(itr))
    end
end

function iterate(c::Combinations,state::Int64=0)
    if state>=length(c)
        return nothing
    end
    indices=digits(state,base=c.itrsize,pad=c.count)
    [c.itr[i] for i in (indices .+1)],state+1
end

function length(c::Combinations)
    length(c.itr) ^ c.count
end

###### MAIN

INSTANCE = ARGS[1]

SEED = parse(Int64, ARGS[2])

MAX_IT = parse(Int64, ARGS[3])

Mt = MersenneTwister(SEED)

CONSTRAINED = false

PERC = .5

EDGES_FILE = "data/" * INSTANCE * ".link"
LABEL_FILE = "data/" * INSTANCE * ".label"

OUTPUT_FILE = "out/" * INSTANCE * "-" * string(SEED) * ".txt"

io = open(OUTPUT_FILE, "w")
close(io)

# tolerance epsilon
Tol = 1e-4

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

G = SimpleWeightedGraph(n)

# load dataset
df_edges = CSV.read(EDGES_FILE; header = false)

# matrix of samples attributes
E = convert(Matrix, df_edges)

global M = 0
for j = 1:size(E)[1]
    global M
    M += E[j, 3]
    if E[j, 1] != E[j, 2]
        M += E[j, 3]
    end
    add_edge!(G, E[j, 1], E[j, 2], E[j, 3])
end

# total number of edges
M = .5*M

data = Data(n, 2, k, G)

# SBM constant
C = M*(log(2*M) - 1)

# (k_i k_j)/2M constant matrix for the DC-SBM
# Q = (data.degree * transpose(data.degree))/(2*M)

ground = Solution(data, copy(label))
w_ground = get_omega_DCSBM(ground.m, ground.kappa)

run(MAX_IT)

# perm = collect(Combinations([1,2], n))

# for x in perm
#     sol = Solution(data, x)
#     ll = ll_DCSBM(data, sol.m, sol.kappa)
#     println(x, ll)
# end

# ground_truth = Solution(data, copy(label))
# omega = get_omega_DCSBM(ground_truth.m, ground_truth.kappa)
# [ println("w(", i, "): ", omega[i,:]) for i = 1:data.k ]
# println("ll = ", ll_DCSBM(data, ground_truth.m, ground_truth.kappa))
