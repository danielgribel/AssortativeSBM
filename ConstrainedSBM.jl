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
using Metis
using SparseArrays

import Base.iterate,Base.length

include("Data.jl")

mutable struct Solution
    # Log-likelihood value
    ll::Float64

    # Dataset
    data::Data
    
    # Sample-community assignment -- class representation
    y::Array{Int}
    
    # Sample-community assignment -- binary indicator
    z::Array{Int}

    # Number of edges for a pair of clusters
    m::Array{Int}

    # Sum of sample degrees on each cluster
    kappa::Array{Int}
end

function Solution(data, y)
    # Initialize cluster binary-indicator variable
    z = zeros(Int, data.n, data.k)
    [ z[i, y[i]] = 1 for i = 1:data.n ]
    
    # Initialize the number of edges for a pair of clusters
    m = compute_m(data, y)
    
    # Initialize the sum of degrees of clusters
    kappa = compute_kappa(data, y)
    
    # Calculate SBM probability matrix
    omega = get_omega_DCSBM(m, kappa)

    # Check assortativity conditions and calculate the log-likelihood
    if is_assortative(omega)
        ll = ll_DCSBM(data, m, kappa)
    else
        ll = solve_DCSBM(data, m, kappa, omega)
    end

    # Create solution
    return Solution(ll, data, y, z, m, kappa)
end

# Update matrix of the number of edges going from cluster r to s
function compute_m(data, y)
    G = data.G
    m = zeros(Int, data.k, data.k)
    for i = 1:data.n
        m[ y[i], y[i] ] += Int(G.weights[i, i])
        for j = (i+1):data.n
            m[ y[i], y[j] ] += Int(G.weights[i, j])
            m[ y[j], y[i] ] += Int(G.weights[j, i])
        end
    end
    return m
end

# Update array of sum of degrees in each cluster
function compute_kappa(data, y)
    kappa = zeros(Int, data.k)
    [ kappa[ y[i] ] += data.degree[i] for i = 1:data.n ]
    return kappa
end

# SBM log-likelihood computation
function ll_SBM(data, m, q)
    k = data.k
    ll = 0.
    [ ll -= .5 * (m[r, r] > 0) * m[r, r] * log(max(1e-10, m[r, r]/(q[r] * q[r]))) for r = 1:k ]
    [ ll -= (m[r, s] > 0) * m[r, s] * log(max(1e-10, m[r, s]/(q[r] * q[s]))) for r = 1:k for s = (r+1):k ]
    return ll - C
end

# DC-SBM log-likelihood computation
function ll_DCSBM(data, m, kap)
    k = data.k
    ll = 0.
    [ ll -= .5 * (m[r, r] > 0) * m[r, r] * log(max(1e-10, m[r, r]/(kap[r] * kap[r]))) for r = 1:k ]
    [ ll -= (m[r, s] > 0) * m[r, s] * log(max(1e-10, m[r, s]/(kap[r] * kap[s]))) for r = 1:k for s = (r+1):k ]
    return ll - C
end

function solve_SBM(data, m, q, omega)
    k = data.k
    
    # SBM parameters
    w = Variable(k, k)
    w.value = omega

    # Threshold variable
    lambda = Variable()

    # Objective function
    f  = 0.5*sum( ( m[r, r]*log(w[r, r]) - q[r]*q[r]*w[r, r] ) for r = 1:k )
    f += sum( ( m[r, s]*log(w[r, s]) - q[r]*q[s]*w[r, s] ) for r = 1:k for s = (r+1):k )

    problem = maximize(f)

    # Add feasibility constraints: w[r, s] > 0 and w[r, s] == w[s, r]
    add_feasibility_constraints(problem, w, k)

    # Add model constraints
    add_strong_constraints(problem, w, lambda, k)

    solve!(problem, ECOS.Optimizer(verbose = false))

    # println("Optval = ", -problem.optval)
    # println("Omega = ", w.value)

    return -problem.optval
end

function solve_DCSBM(data, m, kappa_, omega)
    k = data.k
    
    T = zeros(Float64, k, k)

    # Degrees term
    [ T[r,s] = (kappa_[r]*kappa_[s])/(2*M) for r = 1:k for s = r:k ]
    [ T[s,r] = T[r,s] for r = 1:k for s = (r+1):k ]
    
    # SBM parameters
    w = Variable(k, k)
    w.value = omega

    # Threshold variable
    lambda = Variable()

    # Objective function
    f  = 0.5*sum( ( m[r,r]*log(w[r,r]) - w[r,r]*T[r,r] ) for r = 1:k )
    f += sum( ( m[r,s]*log(w[r,s]) - w[r,s]*T[r,s] ) for r = 1:k for s = (r+1):k )

    problem = maximize(f)

    # Add feasibility constraints: w[r,s] > 0 and w[r,s] == w[s,r]
    add_feasibility_constraints(problem, w, k)

    # Add model constraints
    add_strong_constraints(problem, w, lambda, k)

    # Solve the constrained convex problem
    solve!(problem, ECOS.Optimizer(verbose = false))
    
    # println("Optval = ", -problem.optval)
    # println("Omega = ", w.value)
    
    return -problem.optval
end

# Add constraints related to feasibility
function add_feasibility_constraints(problem, w, k)
    [ problem.constraints += [ w[r, s] >= Tol ] for r = 1:k for s = r:k ]
    [ problem.constraints += [ w[r, s] == w[s, r] ] for r = 1:k for s = (r+1):k ]
end

# Add strong assortative constraints
function add_strong_constraints(problem, w, lambda, k)
    [ problem.constraints += [ w[r, r] >= (lambda + Tol) ] for r = 1:k ]
    [ problem.constraints += [ w[r, s] <= (lambda - Tol) ] for r = 1:k for s = (r+1):k ]
end

# Add weak assortative constraints
function add_weak_constraints(problem, w, k)
    [ problem.constraints += [ w[r, r] >= (w[r, s] + Tol) ] for r = 1:k for s = (r+1):k ]
    [ problem.constraints += [ w[r, r] >= (w[s, r] + Tol) ] for r = 1:k for s = (r+1):k ]
end

# Check strong assortativity condition
function is_strongly_assortative(w)
    return minimum(diag(w)) > maximum(w - Diagonal(w))
end

# Check weak assortativity condition
function is_weakly_assortative(w)
    k = size(w)[1]
    for r = 1:k
        for s = 1:k
            if r != s && (w[r,r] <= w[r,s])
                return false
            end
        end
    end
    return true
end

# Check relaxed assortativity (c-assortative) condition
function is_relaxed_assortative(w, PERC)
    k = size(w)[1]
    count = 0
    for r = 1:k
        arr = findall(w[r, :] .== maximum(w[r, :]))
        if length(arr) == 1 && r == arr[1]
            count += 1
        end
    end
    return (count/k) >= PERC
end

# Check assortativity conditions
function is_assortative(w)
    if CONSTRAINED
        return is_strongly_assortative(w)
    end
    return true
end

function eval_relocate(sol, p, tgt)
    # source cluster
    src = sol.y[p]

    m_ = copy(sol.m)

    for v in neighbors(sol.data.G, p)
        e = Int(sol.data.G.weights[p, v])
        if p == v
            m_[src, src] -= e
            m_[tgt, tgt] += e
        else
            m_[sol.y[v], src] -= e
            m_[src, sol.y[v]] -= e
            m_[tgt, sol.y[v]] += e
            m_[sol.y[v], tgt] += e
        end
    end

    kappa_ = copy(sol.kappa)
    kappa_[ src ] -= sol.data.degree[p]
    kappa_[ tgt ] += sol.data.degree[p]

    sbm_cost = ll_DCSBM(sol.data, m_, kappa_)

    if sbm_cost < sol.ll
        omega = get_omega_DCSBM(m_, kappa_)
        if is_assortative(omega)
            update_ll(sol, sbm_cost)
            update_param(sol, m_, kappa_)
            update_assignment(sol, p, src, tgt)
        else
            sbm_cost = solve_DCSBM(sol.data, m_, kappa_, omega)
            if sbm_cost < sol.ll
                update_ll(sol, sbm_cost)
                update_param(sol, m_, kappa_)
                update_assignment(sol, p, src, tgt)
            end
        end
    end
end

# Update the log-likelihood value
function update_ll(sol, ll)
    sol.ll = ll
end

# Update the SBM parameters
function update_param(sol, m, kappa)
    sol.m = m
    sol.kappa = kappa
end

# Update the solution assignment
function update_assignment(sol, i, src, tgt)
    sol.z[i, src] = 0
    sol.z[i, tgt] = 1
    sol.y[i] = tgt
end

# Calculate maximum-likelihood of \omega in the SBM model, given z
function get_omega_SBM(m, q)
    k = length(q)
    w = [ (m[r, s])/(q[r] * q[s]) for r = 1:k, s = 1:k ]
    return w
end

# Calculate maximum-likelihood of \omega in the DC-SBM model, given z
function get_omega_DCSBM(m, kappa)
    k = length(kappa)
    w = [ 2.0*M*(m[r, s])/(kappa[r] * kappa[s]) for r = 1:k, s = 1:k ]
    return w
end

function localsearch(sol)
    prev_ll = Inf
    curr_ll = sol.ll
    it_changed = true
    
    while ((prev_ll - curr_ll) > Tol) && it_changed
        rdm_samples = randperm(Mt, sol.data.n)
        it_changed = false
        for i in rdm_samples
            prev = sol.y[i]
            rdm_clusters = randperm(Mt, sol.data.k)
            for c in rdm_clusters
                if prev != c
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
    # Random permutation
    rdm_order = randperm(Mt, data.n)

    # Create equaly-sized clusters from the random permutation 
    y = zeros(Int, data.n)
    [ y[ rdm_order[i] ] = ceil(data.k*i/data.n) for i = 1:data.n ]
    
    return y
end

function run(max_it)
    best_ll = Inf
    best_solution = nothing

    for i = 1:max_it
        # Create the initial assignment
        y = initial_assignment(data)

        # Create the initial solution
        sol = Solution(data, y)
        
        # Start CPU time measurement
        t1 = time_ns()
        
        # Apply the local search
        localsearch(sol)
        
        # Finish CPU time measurement
        elapsed_time = (time_ns() - t1)/1.0e9
        
        # Calculate the NMI
        nmi = mutual_information(sol.y, label; normalize = true)
        
        # Calculate the C-rand index
        crand = randindex(sol.y, label)[1]

        # Obtain \omega values
        w = get_omega_DCSBM(sol.m, sol.kappa)

        # Check if the current solution is better than the best solution
        if sol.ll < best_ll
            best_solution = sol
            best_ll = sol.ll
        end

        count_csbm = 0
        for r = 1:data.k
            arr_csbm = findall(w[r, :] .== maximum(w[r, :]))
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
        # [ line *= string(w[r,s]) * " " for r = 1:data.k for s = 1:data.k ]
        line *= "\n"
        print(line)
        # write_output(line)
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
        new{T}(itr, Int64(count), length(itr))
    end
end

function iterate(c::Combinations,state::Int64=0)
    if state >= length(c)
        return nothing
    end
    indices = digits(state, base = c.itrsize, pad = c.count)
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

EDGES_FILE = "data/" * INSTANCE * ".link"
LABEL_FILE = "data/" * INSTANCE * ".label"
OUTPUT_FILE = "out/" * INSTANCE * "-" * string(SEED) * ".txt"

io = open(OUTPUT_FILE, "w")
close(io)

# tolerance epsilon
Tol = 1e-4

label = Int[]

open(LABEL_FILE) do file
    [ push!(label, parse(Int, ln)) for ln in eachline(file) ]
end

# Number of nodes
n = length(label)

# Number of clusters
k = length(unique(collect(label)))

# Create a graph instance
G = SimpleWeightedGraph(n)

# Load the dataset
df_edges = CSV.read(EDGES_FILE; header = false)

# Matrix of edges
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

# Total number of edges
M = .5*M

# Create a data instance
data = Data(n, 2, k, G)

# SBM constant
C = M*(log(2*M) - 1)

CONSTRAINED = true

truth = Solution(data, label)
w_truth = get_omega_DCSBM(truth.m, truth.kappa)

run(MAX_IT)
