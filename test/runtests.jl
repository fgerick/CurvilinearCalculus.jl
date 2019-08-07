using CurvilinearCalculus, LinearAlgebra, StaticArrays, SymPy
using Test



@testset "Spherical coordinates" begin include("spherical.jl") end

@testset "Non-orthogonal elliptical coordinates" begin include("ellipse.jl") end
