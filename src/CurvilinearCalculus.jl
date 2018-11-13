__precompile__(false)

module CurvilinearCalculus

using SymPy, PyCall, LinearAlgebra, StaticArrays #, coordinates

export CoordinateSystem, GenericCoordinates, coordinates, isorthogonal, assume, Q, global_assumptions, dot



abstract type CoordinateSystem end


#convenient definitions for SymPy use:
π = PI
∂ = diff
import LinearAlgebra.dot
symdot(x::Sym,y::Sym) = sum(x.*y)
dot(x::Sym, y::Sym) = symdot(x,y)

@pyimport sympy as sp

global_assumptions=sp.assumptions[:assume][:global_assumptions];

function assumeglobal!(as,global_assumptions=sp.assumptions[:assume][:global_assumptions])
    global_assumptions[:add](as)
end

function assume(var,condition)
    assumeglobal!(eval(:(Q.$(condition)($(var)))),CurvilinearCalculus.global_assumptions)
end
# macro assume(var,condition)
#     :(assumeglobal!(Q.$(condition)($(var)),CurvilinearCalculus.global_assumptions))
# end

const Coordinate = Sym


# const CartesianVector = SVector{3,Sym}

const Vector3D = SVector{3,Sym}
const CoordinateMapping = SVector{3,Sym}

const Metric = SMatrix{3,3,Sym}


macro coordinates(x,y,z)
    x,y,z = eval(:(@syms $x $y $z real=true))
    q = SVector(x,y,z)
    return q
end

struct GenericCoordinates <: CoordinateSystem
    q::SVector{3,Coordinate}
    g_cov::Vector{Vector3D}
    g_contra::Vector{Vector3D}
    G::Metric
    invG::Metric
    J::Sym
    Γ::SArray{Tuple{3,3,3},Sym}

    function GenericCoordinates(r::CoordinateMapping, q::SVector{3,Coordinate})
        g = [∂.(r,q) for q in q]
        J = refine(g[1] ⋅ (g[2] × g[3]))
        g¹= simplify.(1/J * g[2] × g[3])
        g²= simplify.(1/J * g[3] × g[1])
        g³= simplify.(1/J * g[1] × g[2])
        g_contra = [g¹, g², g³]
        G = simplify.([g[i] ⋅ g[j] for i = 1:3, j = 1:3])
        # invG = simplify.(inv(G))
        invG = simplify.([g_contra[i] ⋅ g_contra[j] for i = 1:3, j = 1:3])
        Γ =  simplify.([sum([invG[i,p]*(∂(G[p,j],q[k]) + ∂(G[p,k],q[j]) - ∂(G[j,k],q[p])) for p=1:3])/2 for i=1:3,j=1:3,k=1:3]);
        return new(q,g,g_contra,G,invG,J,Γ)
    end
end

isorthogonal(C::GenericCoordinates) = isdiag(C.G)

struct CovariantVector
    r::Vector3D
    G::GenericCoordinates
end

struct ContravariantVector
    r::Vector3D
    G::GenericCoordinates
end

function CovariantVector(x::ContravariantVector)
    return CovariantVector(x.G.invG*x.r,x.G)
end

function ContravariantVector(x::CovariantVector)
    return ContravariantVector(x.G.invG*x.r,x.G)
end

end # module
