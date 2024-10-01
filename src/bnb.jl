"""
    @File: bnb.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""

module BNB

using Polyhedra: func
export solveBNB, BNBSolver, Logger

##################################################
# Imports and includes, constants and globals
##################################################

using JuMP
using CPLEX
using Test
using Printf
using DataStructures
using LinearAlgebra

using Random
using TimerOutputs
using LinearAlgebra
using Statistics
using CPUTime

import Base.+
import Base.-
import Base.*
import Base./

const MOI = JuMP.MathOptInterface
const EPS = 1e-6
const EPS_P = 1e-6
const SHIFT = 1.1
const R = 1e-3
MAX = typemax(Float64)

include("formatted_prints.jl")
include("structures.jl")
include("geometry.jl")
include("model.jl")
include("gtspn_solver.jl")

##################################################
# The main solver class
##################################################
mutable struct BNBSolver{T}
    goals::Vector{T}
    dist_mat::Any # GTSPN matrix

    best::PartialSolution
    lb::PartialSolution

    open

    to

    times_eval::Array{Any}
    lbs_eval::Array{Any}
    ubs_eval::Array{Any}
    nbr_created_nodes::Int64
    nbr_expanded_nodes::Int64
    nbr_cplex_nodes::Int64
    ub_value::Float64
end

mutable struct Logger
    ubs::Vector{Float64}
    lbs::Vector{Float64}
    ctimes::Vector{Float64}

    nbr_created_nodes::Int64
    nbr_expanded_nodes::Int64
    nbr_exact_nodes::Int64

    ubs_sequences::Vector{Int64}
    lbs_sequences::Vector{Int64}

    ts::Vector{Any}
end


function BNBSolver(filename::String, type::String, dim::String)
    targets, n = load(filename, type, dim)
    ub = read_ub_from_file(filename, type, dim)
    solver = BNBSolver(targets, ones(n, n), PartialSolution(Inf), PartialSolution(0.0), 
            nothing, nothing, [], [], [], 0, 0, 0, ub)
    return solver
end 

function Logger()
    return Logger([], [], [], 0, 0, 0, [], [], [])
end

##################################################
# Log functions
##################################################
function log(solver, logger, start_time)
    now_time = (CPUtime_us()  - start_time) * 1e-6 # miliseconds
    push!(logger.ctimes, now_time)
    push!(logger.lbs, solver.lb.reward)
    push!(logger.ubs, solver.ub.reward)

    push!(logger.ubs_sequences, length(solver.ub.sequence))
    push!(logger.lbs_sequences, length(solver.lb.sequence))
end

function log_ts(logger, ts)
    push!(logger.ts, ts)
end

##################################################
function load(filename::String, type::String, dim::String)
    if !isfile(filename)
        return [], 0
    end
   if type == "GTSPN"
        return load_gtspn(filename)
    end
end

function read_ub_from_file(problem::String,type::String, dim::String)
    filename = "ubs-values/$(lowercase(type))_$(lowercase(dim)).log"
    if !isfile(filename)
        return Inf
    end

    for line in readlines(open(filename))
        pairs = split(line,";")
        name = split(pairs[1],":")[2]
        value = parse(Float64, "$(split(pairs[2],":")[2])")
        if occursin(name, problem)
            return value
        end
    end
    Inf
end

@fastmath function tour_length(tour; closed=true)
    total_dist = 0.0
    for i in 1:length(tour) - 1
        total_dist += dist(tour[i], tour[i + 1])
    end
    if closed 
        total_dist += dist(tour[1], tour[length(tour)])
    end
    return total_dist
end

function log(solver, start_time)
    now_time = (CPUtime_us()  - start_time) * 1e-3 #miliseconds
    push!(solver.times_eval, (@sprintf "%0.7f" now_time))
    push!(solver.lbs_eval, solver.lb.length)
    push!(solver.ubs_eval, solver.best.length)
end

##################################################
# BNB methods
##################################################
function select_root(solver, targets::Vector{<:AbstractTarget})
    n = length(targets)
    sequence::Vector{Int64} = [1 for i in 1:3]
    max_len = 0.0
    for i in 1:n
        for j in i + 1:n
            for k in j + 1:n
                len = solver.dist_mat[i,j] + solver.dist_mat[j,k] + solver.dist_mat[k,i]
                if len > max_len 
                    max_len = len
                    sequence = [i, j, k]
                end
            end
        end
    end
    return sequence
end

function compute_bounds(solver, sequence::Vector{Int64}; compute_ub = true)
    node = Node(PartialSolution(0.0, copy(sequence), []), PartialSolution(Inf), false, [], true)
    @timeit solver.to "ComputeLowerBound" node.lb.points, node.lb.length, _  = find_shortest_tour(solver, solver.goals[sequence]; ub=solver.ub_value)
    if length(node.lb.points) == 0
        return Node()
    end
    if compute_ub
        @timeit solver.to "ComputeUpperBound" node.ub = compute_upper_bound(solver, solver.goals, node)
    else
        @timeit solver.to "LoadUpperBound" node.ub = PartialSolution(solver.ub_value,  [], [])
    end
    return node
end

function compute_estimation(solver, sequence::Vector{Int64})
    node = Node(PartialSolution(0.0, copy(sequence), []), PartialSolution(Inf), false, [], false)
    @timeit solver.to "ComputeLowerBoundEstimation" node.lb.points, node.lb.length = estimate_shortest_tour(solver.goals[sequence])
    return node
end

##################################################
# Main solver
##################################################
function solveBNB(filename::String, max_time::Float64, update_fce, dims::String, type::String, config_as_dict)
    # @show filename, max_time, dims, type
    solver = BNBSolver(filename, type, dims)
    logger = Logger()

    solver.to = TimerOutput()
    solver.dist_mat = distance_matrix(solver.goals)

    big_m_values(solver.goals)

    init_time = CPUtime_us() * 1e-3 # [ms]
    start_time = CPUtime_us()
    solver.open = PriorityQueue()

    root::Node = Node()

    @timeit solver.to "Root" begin
        compute_ub = false
        if solver.ub_value == Inf
            compute_ub = true
        end
        @show compute_ub
        @timeit solver.to "SelectSequence" sequence = select_root(solver, solver.goals)
        @timeit solver.to "ComputeBounds" root = compute_bounds(solver, sequence; compute_ub = compute_ub)
        @assert root !== nothing

        solver.lb = root.lb
        solver.best = root.ub
        if solver.best.length < solver.ub_value 
            solver.ub_value = solver.best.length
        end
        update_fce(solver, logger, "Root")
        log(solver, init_time)

    end

    solver.open[root] = root.lb.length   
    solver.nbr_created_nodes += 1

    while (CPUtime_us() * 1e-3 - init_time) * 1e-3 <= max_time
        if length(solver.open) == 0
            @error "Empty queue"
            @assert length(solver.open) == 0
        end
        cur::Node = dequeue!(solver.open)   
        
        if !cur.cplex
            @timeit solver.to "ComputeExact" begin
                @timeit solver.to "ComputeBounds" cur = compute_bounds(solver, cur.lb.sequence) 
                if length(cur.lb.points) == 0                
                    @warn "CPLEX could not find solution"
                    continue
                end
                if cur.ub < solver.best
                    solver.best = cur.ub 
                    if solver.best.length < solver.ub_value 
                        solver.ub_value = solver.best.length
                    end
                    update_fce(solver, logger,"Updated UB - CPLEX")
                    log(solver, init_time)
                end         
                if cur.lb.length <= solver.best.length + EPS
                    # @info "Added new child for branching: $(cur.lb.sequence) (CPLEX lb=$(cur.lb.length))" 
                    solver.open[cur] = cur.lb.length
                    solver.nbr_cplex_nodes += 1
                end
                continue
            end
        end
        solver.nbr_expanded_nodes += 1
        if cur.lb.length + 1e-3 < solver.lb.length
            @error string("Incorrect lower bound ", cur.lb.length, " < ",  solver.lb.length)
            @assert cur.lb.length >= solver.lb.length
        end

        @timeit solver.to "ComputeNotCovered" cur.notCovered = compute_not_covered(solver.goals, cur.lb)    

        solver.lb = cur.lb    
        update_fce(solver, logger,"Updated LB - dequeued")
        log(solver,   init_time)        

        if isempty(cur.notCovered) 
            if cur.ub < solver.best
                solver.best = cur.ub
                update_fce(solver, logger, "Updated UB - Empty notCovered")
                log(solver, init_time)
            end
            update_fce(solver,logger, "Finish - Empty notCovered")
            log(solver, init_time)
            return solver, logger
        else
            @timeit solver.to "Branching" begin 
                t = findmax([x.dist for x in cur.notCovered])[2]
                nnode::Int64 = cur.notCovered[t].label
                
                for i in 1:length(cur.lb.sequence)
                    for n in cur.lb.sequence            
                        @test n != nnode
                    end

                    sequence = copy(cur.lb.sequence)
                    insert!(sequence, i, nnode)
                    @timeit solver.to "ComputeEstimates" child = compute_estimation(solver, sequence)
                    solver.open[child] = child.lb.length
                    solver.nbr_created_nodes += 1
                    # @info "Added new child for branching: $(child.lb.sequence) (lb=$(child.lb.length))" 
                end
            end
        end
        @timeit solver.to "FilterOpen" solver.open = filter!(x->x[1].lb.length <= solver.best.length + EPS, solver.open)
    end
  
    update_fce(solver,  logger, "Final - timeout") 
    log(solver, init_time)
  
    return solver, logger
end

end
