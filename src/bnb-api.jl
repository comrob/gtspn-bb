#!/usr/bin/env julia

"""
    @File: bnb-api.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""

module BNBAPI

################################################
# Imports and includes
################################################

include("bnb.jl")
using .BNB

using ConfParser
using CPUTime
using Printf
using PrettyTables
using DataFrames
using CSV

include("formatted_prints.jl")

df = DataFrame(PROBLEM=String[], METHOD=String[], TYPE=String[], DIM=String[],
    TRIAL=Int64[], LB_LENGTH=Float64[], UB_LENGTH=Float64[], GAP=Float64[],
    CTIME=Float64[])

###############################################
# Config
################################################
config = ConfParse("$(@__DIR__)/bnb.ini")
parse_conf!(config)
verbose = retrieve(config, "info", "verbose") == "true" ? true : false
profiling = retrieve(config, "debug", "profiling") == "true" ? true : false

function parse_to_dict()
    dict = Dict()
    for (_, d) in config._data
        for (key, value) in d
            dict[key] = value[1]
        end
    end
    return dict
end

################################################
# Solve functions
################################################
function solve(filename::String, problem::String, trial::Int64, t_max::Float64, type::String, dimension::Int64, config_as_dict::Any)
    print_info("Solving $(problem) in $(dimension)D, trial no. $(trial), T_max = $(t_max) [s]")

    t_start = CPUtime_us()

    function update_function(solver::BNBSolver, logger::Logger, msg::String)
        t = CPUtime_us() # to seconds
        gap = (solver.best.length - solver.lb.length) /  solver.best.length * 100.0
        s = @sprintf("CPUTime: %6.2f [s] | UB: %6.2f LB: %6.2f GAP: %6.2f Open: %4d - %s", (t - t_start) * 1e-6, solver.best.length, solver.lb.length, gap, length(solver.open), msg)
        if config_as_dict["verbose"] == "true"
            print_update(s)
        end
    end

    solver, logger = solveBNB(filename, t_max, update_function, "$(dimension)D", type, config_as_dict)
  
    t_end = CPUtime_us()
    gap = (solver.best.length - solver.lb.length) / solver.best.length * 100.0

    print_success_message("Problem $(problem) in $(dimension)D trial no. $(trial) -- SOLVED")

    problem_name = "$(problem)"
    if config_as_dict["save-solutions"] == "true"
        if !isdir(config_as_dict["solutions-dir"])
            mkdir(config_as_dict["solutions-dir"])
        end
        if !isdir("$(config_as_dict["solutions-dir"])/$(problem_name)")
            mkdir("$(config_as_dict["solutions-dir"])/$(problem_name)")
        end
        if !isdir("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)")
            mkdir("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)")
        end
        if !isfile("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/configurations.txt")
            touch("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/configurations.txt")
        end

        open("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/configurations.txt", "w") do f
            for c in solver.best.points
                line = ""
                for coord in c.coords
                    line *= "$(coord) "
                end

                write(f, "$(line)\n")
            end
        end

        open("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/sequence.txt", "w") do f
            for c in solver.best.sequence
                write(f, "$(c)\n")
            end
        end

        open("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/configurations-lb.txt", "w") do f
            for c in solver.lb.points
                line = ""
                for coord in c.coords
                    line *= "$(coord) "
                end
                write(f, "$(line)\n")
            end
        end

        open("$(config_as_dict["solutions-dir"])/$(problem_name)/$(trial)/sequence-lb.txt", "w") do f
            for c in solver.lb.sequence
                write(f, "$(c)\n")
            end
        end
    end

    println(solver.to)
    return true, (t_end - t_start) * 1e-6, solver.lb.length, solver.best.length, gap
end


################################################
# Main
################################################
function main()
    filename = retrieve(config, "problem", "problem")
    dir = retrieve(config, "problem", "batch-dir")
    type = retrieve(config, "problem", "type")
    trials = parse(Int64, retrieve(config, "problem", "trials"))
    dimension = parse(Int64, replace(retrieve(config, "problem", "dim"), "D" => ""))
    t_max = parse(Float64, retrieve(config, "problem", "timeout"))

    config_as_dict = parse_to_dict()

    data::Matrix{Float64} = ones(trials, 5)
    if config_as_dict["batch-mode"] == "true"

        # run first for everything to load
        if config_as_dict["first-run"] == "true"
            filename = readdir(dir)[1]
            filename_ = split(filename, "/")[end]
            problem = replace(basename(filename_), ".txt" => "")
            _, time, lb, ub, gap = solve("$(dir)/$(filename)", problem, 1, 3.0, type, dimension, config_as_dict)
        end

        budgets_lengths = Dict()

        for filename in readdir(dir)
            for trial = 1:trials
                filename_ = split(filename, "/")[end]
                problem = replace(basename(filename_), ".txt" => "")
                _, time, lb, ub, gap = solve("$(dir)/$(filename)", problem, trial, t_max, type, dimension, config_as_dict)
                data[trial, :] = [trial time lb ub gap ]
                push!(df, [problem config_as_dict["method"] uppercase(type) "$(dimension)D" trial lb ub gap time])

            end
        end        
    else
        if config_as_dict["rci-run"] == "true"
            print_info("Running on RCI, taking the instance file from command... $(ARGS)")
            filename = ARGS[1]
            
            filename_ = split(filename, "/")[end]
            problem = replace(basename(filename_), ".txt" => "")
            to_solve = "$(dir)/$(filename)"
        else
            filename_ = split(filename, "/")[end]
            problem = replace(basename(filename_), ".txt" => "")
            to_solve = "$(filename)"

        end

      

        if config_as_dict["first-run"] == "true"
            _, time, lb, ub, gap = solve(to_solve, problem, 1, 3.0, type, dimension, config_as_dict)
        end

        for trial = 1:trials
            _, time, lb, ub, gap = solve(to_solve, problem, trial, t_max, type, dimension, config_as_dict)
            data[trial, :] = [trial time lb ub gap]
            push!(df, [problem config_as_dict["method"] uppercase(type) "$(dimension)D" trial lb ub gap time])
        end
    end

    println()
    println(df)

    if retrieve(config, "save", "save-results") == "true"
        if !isdir(retrieve(config, "results", "results-dir"))
            mkdir(retrieve(config, "results", "results-dir"))
        end

        if config_as_dict["rci-run"] == "true"
            filename_ = split(filename, "/")[end]
            problem = replace(basename(filename_), ".txt" => "")
            name = retrieve(config, "results", "results-dir") * "/" * problem * ".csv"
        else
            name = retrieve(config, "results", "results-dir") * "/" * retrieve(config, "results", "results-file")
        end

        CSV.write(name, df, append=true)
    end

end

main()

end
