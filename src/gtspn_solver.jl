"""
    @File: gtspn_solver.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""

##################################################
# GTSPN loader
##################################################

function add(set::TargetSet, neigh::Target)
  neigh.target_label = neigh.label * set.label
  if neigh.has_ellipse && length(neigh.half_planes) > 0
      push!(set.types, HYBRID)
  elseif neigh.has_ellipse && length(neigh.half_planes) == 0
      push!(set.types, ELLIPSE)
  elseif !neigh.has_ellipse
      push!(set.types, POLYHEDRON)
  end
  push!(set.targets, neigh)
end

function loadPolyhedra(lines::Array{String}, idx::Int64, idx_end::Int64)
  start = idx + 1
  max_hp = ((N == 3) ? 12 : 20)
# TODO -- remove ?? 
  half_planes::Vector{HalfPlane} = [HalfPlane() for x in 1:max_hp]
  for i in start:idx_end
      if occursin("b", lines[i])
          start = i
          break
      end
      half_planes[i - start + 1].a = [parse(Float64, x) for x in split(lines[i])]
  end
  for i in start + 1:idx_end
      if length(lines[i]) == 0 
          continue
      end
      half_planes[i - start].b = parse(Float64, split(lines[i])[1])
  end 
  return half_planes
end

function loadEllipse(lines::Array{String}, idx::Int64, idx_end::Int64)
  start = idx
  pm = zeros(N, N) 
  for i in 1:N
      a = split(lines[start + i])
      for j in 1:N
          pm[i,j] = parse(Float64, a[j])
      end
  end

  return pm
end

function loadHybrid(lines::Array{String}, idx::Int64, idx_end::Int64)
  start = idx
  pm = zeros(N, N) 
  for i in 1:N
      a = split(lines[start + i])
      for j in 1:N
          pm[i,j] = parse(Float64, a[j])
      end
  end

  half_planes::Vector{HalfPlane} = [HalfPlane() for x in 1:6]
  start += N + 2
  for i in start:idx_end
      if occursin("b", lines[i])
          start = i
          break
      end
      half_planes[i - start + 1].a = [parse(Float64, x) for x in split(lines[i])]
  end
  for i in start + 1:idx_end
      if length(split(lines[i])) == 0 
          break
      end
      half_planes[i - start].b = parse(Float64, split(lines[i])[1])
  end 

  return half_planes, pm
end

function loadTarget(lines::Array{String}, idx::Int64, idx_end::Int64)
  neigh = Target(parse(Int64, replace(replace(lines[idx + 1], "Q" => ""), ":" => "")))
  start = idx + 2
  type = split(lines[start])[2]

  neigh.qc = Point([parse(Float64, x) for x in split(lines[start + 1])[end - N + 1:end]])
  neigh.ub = Point([parse(Float64, x) for x in split(lines[start + 2])[end - N + 1:end]])
  neigh.lb = Point([parse(Float64, x) for x in split(lines[start + 3])[end - N + 1:end]])

  if type == "Polyhedra"
      neigh.half_planes = loadPolyhedra(lines, start + 4, idx_end - 1)
  elseif type == "Ellipse"
      neigh.pm1 = loadEllipse(lines, start + 4, idx_end - 1)
      neigh.has_ellipse = true
  elseif type == "Hybrid"
      neigh.half_planes, neigh.pm1 = loadHybrid(lines, start + 4, idx_end - 1)
      neigh.has_ellipse = true
  end

  neigh.bounding_box = get_bounding_points(neigh.lb, neigh.ub)
  return neigh
end

function loadTargetSet(lines::Array{String}, idx::Int64, idx_end::Int64, sets::Vector{TargetSet})
  if idx_end < idx return ; end    

  set = TargetSet(parse(Int64, replace(replace(lines[idx + 2], "S" => ""), ":" => "")))
  start = idx + 3

  set.qc = Point([parse(Float64, x) for x in split(lines[start])[end - N + 1:end]])
  D = length(set.qc.coords)
  set.ub = Point([-Inf for x in 1:D])
  set.lb = Point([Inf for x in 1:D])

  neigh_starts = []
  for i in start + 3:idx_end - 1
      if occursin("=", lines[i]) && occursin("Q", lines[i + 1])
          push!(neigh_starts, i)
      end
  end

  if length(neigh_starts) == 1
      a = circshift(neigh_starts, 1)
      add(set, loadTarget(lines, a[1], idx_end - 2))
  else
      for i in length(neigh_starts):-1:1
          a = circshift(neigh_starts, i)
          add(set, loadTarget(lines, a[1], a[2] == neigh_starts[1] ? idx_end - 2 : a[2]))
      end
  end

  for region in set.targets
      for i in 1:D
          if region.lb.coords[i] < set.lb.coords[i]
              set.lb.coords[i] = region.lb.coords[i]
          end
          if region.ub.coords[i] > set.ub.coords[i]
              set.ub.coords[i] = region.ub.coords[i]
          end
      end
  end
  push!(sets, set)

  return sets
end

function load_gtspn(filename::String)
  sets::Vector{TargetSet} = []
  lines = readlines(open(filename))
  set_starts = []
  global N = parse(Int64, split(lines[1])[end])
  for i in 2:(length(lines) - 1)
      if occursin("=", lines[i]) && occursin("=", lines[i + 1])
          push!(set_starts, i)
      end
  end

  for i in length(set_starts):-1:1
      a = circshift(set_starts, i)
      loadTargetSet(lines, a[1], a[2] == set_starts[1] ? length(lines) - 2 : a[2], sets)
  end

  return sets, length(sets)
end

##################################################
# GTSPN functions
##################################################
@inline function distance_matrix(sets::Vector{TargetSet})
    n = length(sets)
    m = length(sets[1].targets)
    distance_matrix = [dist(a.qc, b.qc) for a in sets, b in sets] 
end

function compute_not_covered(sets::Vector{TargetSet}, solution::PartialSolution)
    n = length(sets)
    dists = [NotCovered(i, false, Inf) for i in 1:n]

    for set in sets
        set_covered = false
        if set.label in solution.sequence
            dists[set.label] = NotCovered(set.label, true, 0.0)
            continue
        end
        for region in set.targets
            if !set_covered
                for i in 1:length(solution.points)
                    aperm = circshift(solution.points, i)
                    a = aperm[1]; b = aperm[2];

                    td, pt = dist_point_to_segment(a, b, region.qc)
                    is_covered, dist = check_point_in_region_with_dist(region, pt)
                    if is_covered
                        dists[set.label] = NotCovered(set.label, true, dist)
                        set_covered = true
                    else
                        if dist < dists[set.label].dist
                            dists[set.label] = NotCovered(set.label, false, dist)
                        end
                    end
                end
            end
        end
    end   
    return  filter(x->!x.feasib, dists)
end

function compute_upper_bound(solver, sets::Vector{TargetSet}, node::Node)
    n = length(sets)
    ret = node.lb.sequence
    ret_pt = node.lb.points

    updated = true

    while updated
        updated = false
        notCovered = [x for x in 1:n if !(x in ret)]
        if isempty(notCovered)
            break
        end

        best_dist = Inf
        best_idx = -1
        best_pt = nothing
        best_target_label = -1
        m = length(ret)
        for a in 1:m
            b = mod1(a+1, m)
            for i in 1:length(notCovered)
                set = sets[notCovered[i]]

                td = Inf
                pt = nothing
                for region in set.targets
                    _, x = dist_point_to_segment(ret_pt[a], ret_pt[b], region.qc)    
                    dist, fpt = move_point_to_surface(x, region)
    
                    if dist < td
                        td = dist
                        pt = fpt
                        is_ok = check_point_in_region(region, pt) 
                        @assert is_ok 
                    end            
                end

                # td, pt = dist_point_to_segment(ret_pt[a], ret_pt[b], target.center)
                if td < best_dist
                    best_dist = td
                    best_target_label = notCovered[i]
                    
                    best_idx = a
                    best_pt = pt             
                end
            end
        end
        if best_pt !== nothing
            # @show "Insert $(best_target_label) at pos $(best_idx)"
            ret = circshift(ret, -best_idx)
            push!(ret,best_target_label)
            ret_pt = circshift(ret_pt, -best_idx)
            push!(ret_pt, best_pt)
            @assert length(ret) == length(ret_pt)
            updated = true
        end
    end

    @assert validate_solution_feasibility(sets, ret, ret_pt)
    return PartialSolution(tour_length(ret_pt), ret, ret_pt)
end

function check_point_in_region_with_dist(region::Target, pt::Point)
    norm_dist = dist = Inf
    if region.has_ellipse
        norm_dist = normalized_elipse_distance(region, pt)
        if norm_dist > 1 + EPS
            return false, norm_dist
        end
    end
    if length(region.half_planes) > 0
        in_hp, dist = check_point_in_halfplanes(region.half_planes, pt)
        if !in_hp 
            return false, dist
        end
    end
    true, (dist < norm_dist ? dist : norm_dist)
end

function check_point_in_region(region::Target, pt::Point)
    if region.has_ellipse
        dist = normalized_elipse_distance(region, pt)
        if dist > 1 + EPS
            return false
        end
    end
    if length(region.half_planes) > 0
        in_hp, dist = check_point_in_halfplanes(region.half_planes, pt)
        if !in_hp 
            return false
        end
    end
    true
end

function validate_solution_feasibility(sets::Vector{TargetSet}, sequence::Vector{Int64}, solution::Vector{Point})

    for i in 1:length(sequence)
        set = sets[sequence[i]]
        pt = solution[i]
        ok = false
        for region in set.targets
            if check_point_in_region(region, pt)
                ok = true
            end
        end  
        if !ok
            return false
        end  
    end

    true
end


"""
    Centroid-GTSP+ - ral18
"""
function lio(targets, solution, sequence)
    n = length(solution)
    len = tour_length(solution; closed=false)
    max_steps = 10000
    for _=1:max_steps
        for i=n+1:2*n
            prev = solution[mod1(i-1, n)]
            act = solution[mod1(i, n)]
            next = solution[mod1(i+1, n)]

            act_idx = sequence[mod1(i, n)]
            target = targets[act_idx]

            solution[mod1(i, n)] = optimize_location_3_points(target, prev, act, next)
        end
        len1 = tour_length(solution, closed=false)
        if len1 < len 
            len = len1
        else
            break
        end
    end

    return solution
end


"""
    Centroid-GTSP+ - ral18
"""
function optimize_location_3_points(target, line_a, alt_goal, line_b)
    dim = length(target.center.coords)
    steps = ones(dim) .* 0.1
    q_c = target.center
    diff = alt_goal - q_c
    len_2segments = norm1(line_a - alt_goal) + norm1(line_b - alt_goal)

    while any(steps .> 1e-10)        
        for d=1:dim
            step = steps[d]
            diff_2 = copy(diff)
            diff_2[d] += step

            diff_2, k_min = intersection_on_surface(target, diff_2)
            alt_goal2 = q_c + diff_2
            len_2segments_act = norm1(line_a - alt_goal2) + norm1(line_b - alt_goal2)

            if len_2segments_act < len_2segments 
                len_2segments = len_2segments_act
                alt_goal = alt_goal2
                steps[d] *= 2
                diff = copy(diff_2)
            else
                steps[d]  *=  abs(step) > 1e-10 ? -0.1 : -1
            end
        end
    end
    return alt_goal

end 

##################################################
# GTSPN -Big M methods
##################################################

@inline function calculate_m(region::Target, set::TargetSet) 
    M::Vector{Float64} = []
  
    if region.has_ellipse
        best_vi = 0.0
        rl = 0; rpt = Point();
        for region1 in set.targets
            if region1.label == region.label 
                continue
            end
            for pt in region1.bounding_box
                vi = transpose(pt.coords - region.qc.coords) * region.pm1 * (pt.coords - region.qc.coords)
                if vi > best_vi
                    best_vi = vi
                    rl = region1.label
                    rpt = pt
                end
            end
        end
        region.m = best_vi * SHIFT
    end
    if length(region.half_planes) > 0
        for hp in region.half_planes
            best_vi = 0.0
            rl = 0; rpt = Point();
            for region1 in set.targets
                if region1.label == region.label 
                    continue
                end
                for pt in region1.bounding_box
                    vi = abs(transpose(hp.a) * pt.coords - hp.b)
                    if vi > best_vi 
                        best_vi = vi
                        rl = region1.label
                        rpt = pt
                    end
                end
            end
            hp.m = best_vi * SHIFT
        end
    end
end

function big_m_values(sets::Vector{TargetSet})
    n = length(sets); m = length(sets[1].targets)
    M::Matrix{Array{Float64}} = [[] for i in 1:n, j in 1:m]
    for set in sets
        for region in set.targets
            calculate_m(region, set)
        end  
    end
end

