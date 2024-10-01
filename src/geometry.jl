"""
    @File: geometry.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""

# http://geomalgorithms.com/a02-_lines.html
@inline function dist_point_to_segment(a::Point, b::Point, p::Point) # |ab| = segment
    # a ~ p0; b ~ p1
      v = b - a # p1 - p0
      w = p - a # p - p0
  
      c1 = dot1(w, v)
      if c1 <= 0.0
          td = dist(p, a)
          return td, a
      end
  
      c2 = dot1(v, v)
      if c2 <= c1
          td = dist(p, b)
          return td, b
      end
      pt = a + (c1 / c2) * v
      return dist(p, pt), pt
  end
  
  
  @inline function point_on_plane(a::Array{Float64}, b::Float64, pt::Point; eps::Float64=EPS_P)
      x = dot(a, pt.coords) - b
      return x <= eps, x
  end
  
  @inline function check_point_in_halfplanes(half_planes::Array{HalfPlane}, pt::Point; eps::Float64=EPS_P)
      ok = true
      min_td = Inf
      for hp in half_planes
          d, td = point_on_plane(hp.a, hp.b, pt, eps=eps)
          if !d
              ok = false
              if min_td > td
                  min_td = td
              end
          end
      
      end
      return ok, min_td
  end
  
  @inline function normalized_elipse_distance(el::Target, x::Point)
      diff = x.coords - el.qc.coords
      return transpose(diff) * el.pm1 * diff
  end
  
  @inline function normalized_elipse_distance_diff(el::Target, x::Point)
      diff = x.coords
      return transpose(diff) * el.pm1 * diff
  end
  
  function intersection_on_surface(target::Target, diff_original::Point)
      k_new = Inf
      q_c = target.qc;
  
    # // a . [x,y,z] - b <= 0
    # // a . (q_c + k*diff) - b = 0
    # // (a . q_c) + k*(a . diff) - b = 0
    # // k = (b - (a . q_c)) / (a . diff)
      k_new = Inf 
      for plane in target.half_planes
          a_dot_diff = dot(plane.a, diff_original.coords) 
        
          if a_dot_diff > 1e-5
              k = (plane.b - dot(plane.a, q_c.coords)) / a_dot_diff
              if k < k_new
                  k_new = k;
              end
          end
      end
  
      norm_dst = normalized_elipse_distance_diff(target, diff_original); 
      el_k = 1 / sqrt(norm_dst);
  
      if (el_k < k_new)
          k_new = el_k;
      end
    
      diff_new = diff_original * k_new;
      return diff_new, k_new
  end
  
  function move_point_to_surface(x_coords::Point, region::Target)
      diff = x_coords - region.qc
      diff_new, k_min = intersection_on_surface(region, diff)
      fpt = diff_new + region.qc
  
      if k_min > 1.0
          fpt = x_coords
      end
      @assert check_point_in_region(region, fpt) 
      td = dist(fpt, x_coords)
      return td, fpt
  end