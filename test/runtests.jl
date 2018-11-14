using CurvilinearCalculus, LinearAlgebra, StaticArrays, SymPy
using Test



@time @testset "Spherical coordinates" begin include("spherical.jl") end
