"""
    @File: model.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""


function ellipse_constraint(model, region, p, y; bin=nothing)
    A = convert(Array{Float64,2}, cholesky(region.pm1).U)
    ellipse = transpose(y[:]) * y[:]
    tmp = (p[:] - region.qc.coords)
    tmp2 = A * tmp
    @constraint(model, tmp2 .== y[:])
    if bin === nothing
        @constraint(model, ellipse <= 1.0)
    else
        @constraint(model, ellipse <= 1 + region.m * bin + R)
    end
end

function plane_constraint(model, hp, p; bin=nothing)
    if bin === nothing
        @constraint(model, transpose(hp.a) * p[:] - hp.b <= 0.0)
    else
        @constraint(model, transpose(hp.a) * p[:] - hp.b <= hp.m * bin + R)
    end
end

function plane_constraints(model, region, p; bin=nothing)
    for hp in region.half_planes
        plane_constraint(model, hp, p; bin=bin)
    end
end

function generate_constraints(model, point, set::TargetSet)
    D = length(point)
    m = length(set.targets)

    bin = @variable(model, [1:m], Bin)

    @constraint(model, sum(bin[:]) == m - 1)

    for region in set.targets
        if region.has_ellipse
            y_tmp = @variable(model, [1:D])
            ellipse_constraint(model, region, point, y_tmp; bin=bin[region.label])
        end
        plane_constraints(model, region, point; bin=bin[region.label])
    end
end

function generate_constraints(model, point, region::Target)
    D = length(point)

    @variable(model, y[1:D])

    if region.has_ellipse
        ellipse_constraint(model, region, point, y)
    else
        plane_constraints(model, region, point)
    end
end

function generate_constraints_box(model, point, set::TargetSet)
    D = length(point)
    @constraint(model, point .<= set.ub.coords)
    @constraint(model, point .>= set.lb.coords)
end

function generate_constraints_box(model, point, target::TargetRegion)
    D = length(point)
    @constraint(model, point .<= target.ub.coords)
    @constraint(model, point .>= target.lb.coords)
end

function estimate_shortest_tour(targets::Vector{<:AbstractTarget}; closed::Bool=true)
    n = length(targets)

    # CPX_PARAM_THREADS - no. parallel threads invoked by parallel CPLEX parallel
    # CPX_PARAM_PARALLELMODE - parallel optimization mode (-1 - opportunistic, 0 - Automatic, 1 - deterministic)
    # CPX_PARAM_BARQCPEPCOMP - Convergence tolerance for quadratically constrained problems
    # CPX_PARAM_SCRIND - messages
    # CPX_PARAM_PREIND - presolve
    # CPX_PARAM_REDUCE - Primal and dual reduction  (0 - no reduction, 1 - only primal, 2 - only dual, 3 - both/default)
    # CPX_PARAM_EPRHS - feasibility tolerance
    # CPX_PARAM_TILIM - maximum time, in seconds
    # CPX_PARAM_NUMERICALEMPHASIS - mphasizes precision in numerically unstable or difficult problems.  (0/1)

    model = Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_THREADS=1, CPX_PARAM_PARALLELMODE=1,
        CPX_PARAM_BARQCPEPCOMP=1e-12, CPX_PARAM_SCRIND=0,
        CPX_PARAM_EPRHS=1e-6, CPX_PARAM_NUMERICALEMPHASIS=1
    )
    )  #1e+1 -> spojim stredy, 1e-1 -> body inside of disk
    # model = Model(Ipopt.Optimizer) # ,"tol"<=1e-6,  "print_level"==0, "acceptable_constr_viol_tol"<=1e-3))
    # set_optimizer_attribute(model, "print_level", 0)

    # f -- objective z  
    if closed
        @variable(model, 0 <= f[1:n] <= MAX)
    else
        @variable(model, 0 <= f[1:n-1] <= MAX)
    end

    # dimensionality
    if typeof(targets[1]) <: TargetRegion
        D = length(targets[1].center.coords)
    else
        D = length(targets[1].qc.coords)
    end

    # x, y, z variables - final locations
    @variable(model, x[1:n, 1:D])
    @variable(model, y[1:n, 1:D])

    # temporary for distance calculation of the consecutive points
    @variable(model, w[1:n, 1:D])

    # objective min sum(f_i)
    @objective(model, Min, sum(f))

    for i in 1:length(targets)
        target = targets[i]
        point = x[i, :]
        generate_constraints_box(model, point, target)
    end

    for i = 1:(closed ? n : n - 1)
        @constraint(model, f[i]^2 >= dot(w[i, :], w[i, :]))
    end

    for i = 1:(closed ? n : n - 1)
        @constraint(model, w[i, :] .== x[i, :] - x[mod1(i + 1, n), :])
    end

    try
        JuMP.optimize!(model)
    catch e
        @warn length(targets)
        @warn string("CPLEX Error: ", e)
    end

    if termination_status(model) == MOI.OPTIMAL #  || termination_status(model) == MOI.LOCALLY_SOLVED
        x_coords = JuMP.value.(x)
        if JuMP.has_lower_bound(x[1])
            @warn JuMP.lower_bound(x[1])
        end
        obj = JuMP.objective_value(model)
        points = [Point(x_coords[i, 1:D]) for i in 1:n]

        return points, obj
    end
    @assert false
    return [], 0.0
end


function find_shortest_tour(solver, targets::Vector{<:AbstractTarget}; closed::Bool=true, ub::Float64=100000000.0)
    n = length(targets)

    # CPX_PARAM_THREADS - no. parallel threads invoked by parallel CPLEX parallel
    # CPX_PARAM_PARALLELMODE - parallel optimization mode (-1 - opportunistic, 0 - Automatic, 1 - deterministic)
    # CPX_PARAM_BARQCPEPCOMP - Convergence tolerance for quadratically constrained problems
    # CPX_PARAM_SCRIND - messages
    # CPX_PARAM_PREIND - presolve
    # CPX_PARAM_REDUCE - Primal and dual reduction  (0 - no reduction, 1 - only primal, 2 - only dual, 3 - both/default)
    # CPX_PARAM_EPRHS - feasibility tolerance
    # CPX_PARAM_TILIM - maximum time, in seconds
    # CPX_PARAM_NUMERICALEMPHASIS -emphasizes precision in numerically unstable or difficult problems.  (0/1)
    # CPX_PARAM_OBJULIM = Upper objective value limit

    global model = Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_THREADS=1, CPX_PARAM_PARALLELMODE=1,
        CPX_PARAM_BARQCPEPCOMP=1e-12, CPX_PARAM_SCRIND=0,
        CPX_PARAM_EPRHS=1e-6, CPX_PARAM_NUMERICALEMPHASIS=1,
        CPX_PARAM_OBJULIM=ub, CPXPARAM_TimeLimit=1800
    ))

    # f -- objective z  
    if closed
        @variable(model, 0 <= f[1:n] <= MAX)
    else
        @variable(model, 0 <= f[1:n-1] <= MAX)
    end

    # dimensionality
    if typeof(targets[1]) <: TargetRegion
        D = length(targets[1].center.coords)
    else
        D = length(targets[1].qc.coords)
    end

    # x, y, z variables - final locations
    @variable(model, x[1:n, 1:D])
    @variable(model, y[1:n, 1:D])

    # temporary for distance calculation of the consecutive points
    @variable(model, w[1:n, 1:D])

    # objective min sum(f_i)
    @objective(model, Min, sum(f))

    for i in 1:length(targets)
        target = targets[i]
        point = x[i, :]
        generate_constraints(model, point, target)
    end

    for i = 1:(closed ? n : n - 1)
        @constraint(model, f[i]^2 >= dot(x[i, :] - x[mod1(i + 1, n), :], x[i, :] - x[mod1(i + 1, n), :]))
    end

    try
        JuMP.optimize!(model)
    catch e
        @warn length(targets)
        @warn string("CPLEX Error: ", e)
    end

    x_coords = JuMP.value.(x)
    if JuMP.has_lower_bound(x[1])
        @warn JuMP.lower_bound(x[1])
    end
    obj = JuMP.objective_value(model)
    points = [Point(x_coords[i, 1:D]) for i in 1:n]

    if typeof(targets[1]) <: TargetRegion
        return points, obj, []
    end

    # Check GTSPN validity
    orig_points = points
    for i in 1:length(targets)
        set = targets[i]
        pt = points[i]
        ok = false
        for s in set.targets
            if check_point_in_region(s, pt)
                ok = true
            end
        end

        if !ok
            @timeit solver.to "CPLEXSolutionNotPrecise" begin
                ok2 = false
                new_pt = Point()
                min_dist = Inf
                for region in set.targets
                    td, fpt = move_point_to_surface(pt, region)
                    if td < min_dist
                        min_dist = td
                        new_pt = fpt
                    end
                end
                for s in set.targets
                    if check_point_in_region(s, new_pt)
                        ok2 = true
                    end
                end
                if ok2
                    points[i] = Point(new_pt.coords, true)
                    ok = true
                end
            end
        end
        @assert ok
    end

    return points, tour_length(points), orig_points

end


@inline function get_bounding_points(lb, ub)
    mask = vcat([1 for i in 1:N], [0 for i in 1:N])
    table = []
    for i in 1:2^N
        c = circshift(mask, i)
        push!(table, c[1:N])
    end

    pts = zeros(2^N, N)
    for i in 1:2^N
        for j in 1:N
            pts[i, j] = table[i][j] == 0 ? lb.coords[j] : ub.coords[j]
        end
    end

    return [Point(pts[i, :]) for i in 1:2^N]
end

function get_bound(target::AbstractTarget, goal::MOI.OptimizationSense)
    if typeof(target) <: TargetRegion
        D = length(target.center.coords)
    else
        D = length(target.qc.coords)
    end
    bound = [-Inf for x in 1:D]

    for i in 1:D
        model = Model(with_optimizer(CPLEX.Optimizer, #, CPX_PARAM_THREADS = 1, 
            CPX_PARAM_PARALLELMODE=1, CPX_PARAM_BARQCPEPCOMP=1e-12, CPX_PARAM_SCRIND=0,
            CPX_PARAM_EPRHS=1e-09, CPX_PARAM_PREIND=1)) # CPX_PARAM_SCRIND = 0 

        @variable(model, x[1:D])

        generate_constraints(model, x, target)

        @objective(model, goal, x[i])
        try
            JuMP.optimize!(model)
        catch e
            @warn length(targets)
            @warn string("CPLEX Error: ", e)
        end
        if termination_status(model) == MOI.OPTIMAL  #|| termination_status(model) == MOI.LOCALLY_SOLVED
            if value.(x)[i] > bound[i]
                bound[i] = value.(x)[i]
            end
            if typeof(target) <: Target
                pt = Point(value.(x))
                if target.has_ellipse && length(target.half_planes) > 0
                    continue
                end
                @assert check_point_in_region(target, pt)
            end
        else
            @assert false
        end
    end

    return Point(bound)
end

function lower_upper_bounds(target::AbstractTarget)
    return get_bound(target, MOI.MAX_SENSE), get_bound(target, MOI.MIN_SENSE)
end