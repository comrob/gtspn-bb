"""
    @File: structures.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""
##################################################
# Structures
##################################################
struct Point
  coords::Vector{Float64}
  moved::Bool
end

Point() = Point([], false)
Point(x::Float64, y::Float64) = Point([x, y], false)
Point(x::Float64, y::Float64, z::Float64) = Point([x, y, z], false)
Point(x::Array{Float64}) = Point(copy(x), false)

@inline (+)(a::Point, b::Point) = Point(a.coords .+ b.coords)
@inline (-)(a::Point, b::Point) = Point(a.coords .- b.coords)
@inline (*)(a::Point, b::Real) = Point(a.coords .* b)
@inline (*)(a::Real, b::Point) = (*)(b, a)
@inline (/)(a::Point, b::Real) = Point(a.coords ./ b)

@inline (dot1)(a::Point, b::Point) = dot(a.coords, b.coords)
@inline (norm1)(a::Point) = norm(a.coords)
@inline (angle)(a::Point, b::Point) = acos( dot1(a,b) / (norm1(a)*norm1(b))  )

@inline function dist(a::Point, b::Point)
  diff = a - b
  sqrt(dot1(diff, diff))
end

@inline function dist(a::Point)
  sqrt(dot1(a, a))
end

@inline function dist(a::Array{Float64}, b::Array{Float64})
  diff = a - b
  sqrt(dot(diff, diff))
end

abstract type AbstractTarget end

##################################################
# Structures - CETSP/ETSPN
##################################################
mutable struct TargetRegion <: AbstractTarget # Target region given by the location and its sensing radius
  label::Int64
  center::Point 
  radius::Float64
  lb::Point
  ub::Point
end  

##################################################
# Structures - GTSPN
##################################################
mutable struct HalfPlane
  a::Array{Float64}
  b::Float64
  m::Float64
end

HalfPlane() = HalfPlane([], 0.0, 0.0)

mutable struct Vertex
  hp::Array{Int64}
  f::Int64
end

mutable struct Target <: AbstractTarget
  label::Int64
  qc::Point
  ub::Point
  lb::Point
  has_ellipse::Bool
  pm1::Matrix{Float64} # Array{Float64}
  half_planes::Vector{HalfPlane}
  vertices::Vector{Vertex}
  target_label::Int64
  m::Float64
  bounding_box::Array{Point}
end

Target(label) = Target(label, Point(), Point(), Point(), false, zeros(N, N), [], [], 0, 0.0, [])

@enum TType ELLIPSE POLYHEDRON HYBRID 

mutable struct TargetSet <: AbstractTarget
  label::Int64
  targets::Vector{Target}
  types::Vector{TType}
  ub::Point
  lb::Point
  qc::Point
end

TargetSet(label) = TargetSet(label, [], [], Point(), Point(), Point())

##################################################
# Structures - Branch and Bound
##################################################

mutable struct PartialSolution
  length::Float64
  sequence::Vector{Int64}
  points::Vector{Point}
end

PartialSolution() = PartialSolution(0.0, [], [])
PartialSolution(x::Float64) = PartialSolution(x,[], [])

Base.:<(a::PartialSolution, b::PartialSolution) = a.length < b.length

struct DM
  distances::Matrix{Float64}
end

struct NotCovered
  label::Int64
  feasib::Bool
  dist::Float64
  edge::Vector{Int64}
end

NotCovered(l, f, d) = NotCovered(l, f, d, [])

Base.:<(a::NotCovered, b::NotCovered) = a.dist < b.dist
Base.:isless(a::NotCovered, b::NotCovered) = a.dist < b.dist

mutable struct Node
  lb::PartialSolution
  ub::PartialSolution
  isFeasible::Bool
  notCovered::Vector{NotCovered}
  cplex::Bool
end

Node() = Node(PartialSolution(0.0),PartialSolution(Inf), false, [], false)
# Node(s, f, b, c) = Node(s, PartialSolution(0.0), f, b, c)

mutable struct Branching
  candidates::Vector{Node} 
end