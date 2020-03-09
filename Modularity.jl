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
    ll = ll_modularity(data, z)
    return Solution(ll, data, y, z, m, kappa)
end

function compute_m(data, y)
    G = data.G
    # k = length(unique(collect(y)))
    m = zeros(Int, data.k, data.k)
    for i=1:data.n
        m[ y[i], y[i] ] += Int(G.weights[i, i])
        for j=(i+1):data.n
            m[ y[i], y[j] ] += Int(G.weights[i, j])
            m[ y[j], y[i] ] += Int(G.weights[j, i])
        end
    end
    return m
end

function compute_kappa(data, y)
    # k = length(unique(collect(y)))
    kappa = zeros(Int, data.k)
    [ kappa[ y[i] ] += data.degree[i] for i=1:data.n ]
    return kappa
end

function get_omega_DCSBM(data, y)
    m = compute_m(data, y)
    kappa = compute_kappa(data, y)
    # k = length(kappa)
    k = length(unique(collect(y)))

    B = []
    for i=1:length(kappa)
        if kappa[i] != 0
            push!(B, i)
        end
    end

    m = m[B, B]
    kappa = kappa[B]

    w = [ 2.0*M*(m[r, s])/(kappa[r] * kappa[s]) for r=1:k, s=1:k ]
    return w
end

function ll_modularity(data, z)
    n = data.n
    ll = sum( (T[i, j] *  transpose(z[i,:]) * z[j,:]) for i=1:n for j=1:n )
    ll = ll/(2*M)
    return -ll
end

function eval_relocate(sol, i, t)
    src = sol.y[i]
    z_ = copy(sol.z)
    z_[i, src] = 0
    z_[i, t] = 1
    
    prev_cost = 2.0*sum( (T[i, j] * transpose(sol.z[i,:]) * sol.z[j,:]) for j=1:sol.data.n )
    post_cost = 2.0*sum( (T[i, j] * transpose(z_[i,:]) * z_[j,:]) for j=1:sol.data.n )
    
    prev_cost = -prev_cost/(2*M)
    post_cost = -post_cost/(2*M)

    if post_cost < prev_cost
        # update likelihood
        sol.ll = sol.ll + (post_cost - prev_cost)
        # update assignments
        sol.z[i, src] = 0
        sol.z[i, t] = 1
        sol.y[i] = t
    end
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
    best_ll = Inf
    best_solution = nothing

    for i=1:max_it
        # create initial solution
        y = initial_assignment(data)
        sol = Solution(data, y)
        ll_initial = sol.ll
        t1 = time_ns()
        localsearch(sol)
        elapsed_time = (time_ns() - t1)/1.0e9

        # println(sol.ll, " ", sol.y)

        nmi = mutual_information(sol.y, label; normalize = true)
        crand = randindex(sol.y, label)[1]

        w = get_omega_DCSBM(sol.data, sol.y)
        omega_ii = diag(w)
        omega_ij = w - Diagonal(w)
        
        k = length(unique(collect(sol.y)))

        if sol.ll < best_ll
            best_solution = sol
            best_ll = sol.ll
        end

        count_csbm = 0
        for r = 1:k
            arr_csbm = findall(w[r, :] .== maximum(w[r, :]))
            if length(arr_csbm) == 1 && r == arr_csbm[1]
                count_csbm += 1
            end
        end

        line = INSTANCE * " "
        line *= string(SEED) * " "
        line *= string(k) * " "
        line *= string(sol.ll) * " "
        line *= string(nmi) * " "
        line *= string(count_csbm/k) * " "
        line *= string(elapsed_time) * " "
        for r = 1:k
            for s = 1:k
                line *= string(w[r,s]) * " "
            end
        end
        line *= "\n"
        write_output(line)
    end
    # println(best_solution.y)
end

function write_output(line)
    io = open(OUTPUT_FILE, "a")
    write(io, line)
    close(io)
end

###### MAIN

INSTANCE = ARGS[1]

SEED = parse(Int64, ARGS[2])

MAX_IT = parse(Int64, ARGS[3])

Mt = MersenneTwister(SEED)

CONSTRAINED = false

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

# graph (adjacency matrix)
A = zeros(Int, n, n)

global M = 0
for j = 1:size(E)[1]
    global M
    M += E[j, 3]
    if E[j, 1] != E[j, 2]
        M += E[j, 3]
    end
    add_edge!(G, E[j, 1], E[j, 2], E[j, 3])
    A[E[j, 1], E[j, 2]] = E[j, 3]
    A[E[j, 2], E[j, 1]] = E[j, 3]
end

# total number of edges
M = .5*M

data = Data(n, 2, k, G)

# SBM constant
C = M*(log(2*M) - 1)

# (k_i k_j)/2M constant matrix for the DC-SBM
Q = (data.degree * transpose(data.degree))/(2*M)

T = A - Q

run(MAX_IT)
