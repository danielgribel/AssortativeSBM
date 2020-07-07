using Convex
using ECOS
using Random
using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
using Discreet
using StatsBase
using Clustering
using DelimitedFiles

include("Data.jl")

# Tolerance epsilon
const Tol = 1e-4

global ac_counter
global tt_counter

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
    w = get_omega_DCSBM(m, kappa, data)

    # Check assortativity conditions and calculate the log-likelihood
    if is_assortative(w)
        ll = ll_DCSBM(data, m, kappa)
    else
        ll = solve_convex(data, m, kappa, w)
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

# DC-SBM log-likelihood computation
function ll_DCSBM(data, m, kap)
    k = data.k
    ll = 0.
    [ ll -= .5 * m[r, r] * log(max(1e-10, m[r, r]/(kap[r] * kap[r]))) for r = 1:k ]
    [ ll -= m[r, s] * log(max(1e-10, m[r, s]/(kap[r] * kap[s]))) for r = 1:k for s = (r+1):k ]
    return ll - data.C
end

function solve_convex(data, m, kappa_, omega)
    k = data.k

    T = zeros(Float64, k, k)

    # Degrees term
    [ T[r,s] = (kappa_[r]*kappa_[s])/(2*data.M) for r = 1:k for s = r:k ]
    [ T[s,r] = T[r,s] for r = 1:k for s = (r+1):k ]
    
    # SBM parameters
    w = Variable(k, k)
    w.value = omega

    # Threshold variable
    lambda = Variable()

    # Objective function
    f  = 0.5*sum( ( m[r, r]*log(w[r, r]) - w[r, r]*T[r, r] ) for r = 1:k )
    f += sum( ( m[r, s]*log(w[r, s]) - w[r, s]*T[r, s] ) for r = 1:k for s = (r+1):k )

    # Define a maximization problem
    problem = maximize(f)

    # Add feasibility constraints: w[r,s] > 0 and w[r,s] == w[s,r]
    add_feasibility_constraints(problem, w, k)

    # Add model constraints
    add_strong_constraints(problem, w, lambda, k)

    # Solve the constrained convex problem
    solve!(problem, ECOS.Optimizer(verbose = false))
    
    # println("Optval = ", -problem.optval)
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

# Check assortativity conditions
function is_assortative(w)
    # Constrained-problem flag
    CONSTRAINED = false
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
        omega = get_omega_DCSBM(m_, kappa_, sol.data)
        if is_assortative(omega)
            update_ll(sol, sbm_cost)
            update_param(sol, m_, kappa_)
            update_assignment(sol, p, src, tgt)
        else
            global ac_counter += 1
            sbm_cost = solve_convex(sol.data, m_, kappa_, omega)
            if sbm_cost < sol.ll
                update_ll(sol, sbm_cost)
                update_param(sol, m_, kappa_)
                update_assignment(sol, p, src, tgt)
            end
        end
        global tt_counter += 1
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

# Calculate maximum-likelihood of \omega in the DC-SBM model, given z
function get_omega_DCSBM(m, kappa, data)
    k = length(kappa)
    w = [ 2.0*data.M*(m[r, s])/(kappa[r] * kappa[s]) for r = 1:k, s = 1:k ]
    return w
end

function localsearch(sol, Mt)
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
                # Check sample relocation
                if prev != c
                    eval_relocate(sol, i, c)
                end
                # Update likelihood
                if sol.ll < (curr_ll - Tol)
                    prev_ll = curr_ll
                    curr_ll = sol.ll
                    it_changed = true
                end
            end
        end
    end
end

function initial_assignment(data, Mt)
    # Random permutation
    rdm_order = randperm(Mt, data.n)

    # Create equaly-sized clusters from the random permutation 
    y = zeros(Int, data.n)
    [ y[ rdm_order[i] ] = ceil(data.k*i/data.n) for i = 1:data.n ]
    
    return y
end

function run(max_it, data, label, Mt)
    best_ll = Inf
    best_solution = nothing

    OUTPUT_FILE = "out/" * data.instance * "-" * string(SEED) * ".txt"
    io = open(OUTPUT_FILE, "w")
    close(io)

    for i = 1:max_it
        # Create the initial assignment
        y = initial_assignment(data, Mt)

        # Create the initial solution
        sol = Solution(data, y)

        global ac_counter = 0
        global tt_counter = 0
        
        # Apply the local search
        elapsed_time = @elapsed localsearch(sol, Mt)

        # Calculate the NMI
        nmi = mutual_information(sol.y, label; normalize = true)
        
        # Calculate the C-rand index
        crand = randindex(sol.y, label)[1]

        # Obtain \omega values
        w = get_omega_DCSBM(sol.m, sol.kappa, data)

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

        line = data.instance * " "
        line *= string(SEED) * " "
        line *= string(data.k) * " "
        line *= string(sol.ll) * " "
        line *= string(nmi) * " "
        line *= string(count_csbm/data.k) * " "
        line *= string(elapsed_time) * " "
        line *= string(ac_counter) * " "
        line *= string(tt_counter) * " "
        # [ line *= string(w[r,s]) * " " for r = 1:data.k for s = 1:data.k ]
        line *= "\n"
        print(line)
        write_output(line, OUTPUT_FILE)
    end
end

function write_output(line, OUTPUT_FILE)
    io = open(OUTPUT_FILE, "a")
    write(io, line)
    close(io)
end

###### MAIN
function main(INSTANCE, MAX_IT, SEED)
    Mt = MersenneTwister(SEED)

    EDGES_FILE = "data/" * INSTANCE * ".link"
    LABEL_FILE = "data/" * INSTANCE * ".label"
    
    label = Int[]

    open(LABEL_FILE) do file
        [ push!(label, parse(Int, ln)) for ln in eachline(file) ]
    end

    # Number of nodes
    n = length(label)

    # Number of communities
    k = length(unique(collect(label)))

    # Create a graph instance
    G = SimpleWeightedGraph(n)

    # Load the dataset
    E = readdlm(EDGES_FILE, Int)
    
    # Add edges to the graph
    [ add_edge!(G, E[j, 1], E[j, 2], E[j, 3]) for j = 1:size(E)[1] ]

    # Create a data instance
    data = Data(n, k, G, INSTANCE)
    
    run(MAX_IT, data, label, Mt)
end

INSTANCE = split(ARGS[1], " ")

SEED = parse(Int64, ARGS[2])

MAX_IT = parse(Int64, ARGS[3])

[ main(data, MAX_IT, SEED) for data in INSTANCE ]
