__precompile__(false)

module CurvilinearCalculus

using SymPy, PyCall, LinearAlgebra, StaticArrays, Combinatorics

export CoordinateSystem, GenericCoordinates, coordinates, isorthogonal,
        assume, Q, global_assumptions, dot, Vector3D, CoordinateMapping, Metric,
        CovariantVector, ContravariantVector, CartesianVector

export ∇, grad, div



abstract type CoordinateSystem end

const Coordinate = Sym
const Vector3D = SVector{3,Sym}
const CoordinateMapping = SVector{3,Sym}
const Metric = SMatrix{3,3,Sym}



#Levi-Civita
ϵ = levicivita
ϵt = [ϵ([i,j,k]) for i=1:3,j=1:3,k=1:3]

#convenient definitions for SymPy use:
π = PI
∂ = diff

# define sympy dot as real dot product only to avoid assumptions:
import LinearAlgebra.dot
dot(x::Sym, y::Sym) = sum(x.*y)
dot(x::Vector3D, y::Vector3D) = sum(x.*y)


#make assumptions work (kind of?)
@pyimport sympy as sp
global_assumptions=sp.assumptions[:assume][:global_assumptions];

function assumeglobal!(as,global_assumptions=sp.assumptions[:assume][:global_assumptions])
    global_assumptions[:add](as)
end

function assume(var,condition)
    assumeglobal!(eval(:(Q.$(condition)($(var)))),CurvilinearCalculus.global_assumptions)
end


"""
create coordinate vector
"""
macro coordinates(x,y,z)
    x,y,z = eval(:(@syms $x $y $z real=true))
    q = SVector(x,y,z)
    return q
end

"""
Define Coordinate system...
"""
struct GenericCoordinates <: CoordinateSystem
    q::SVector{3,Coordinate}
    g_cov::Vector{Vector3D}
    g_contra::Vector{Vector3D}
    e_cov::Vector{Vector3D}
    e_contra::Vector{Vector3D}
    G::Metric
    invG::Metric
    J::Sym
    g::Sym
    Γ::SArray{Tuple{3,3,3},Sym}

    function GenericCoordinates(r::CoordinateMapping, q::SVector{3,Coordinate})
        g = [∂.(r,q) for q in q]
        e_cov = [g/norm(g) for g in g]
        J = simplify(g[1] ⋅ (g[2] × g[3]))
        g¹= simplify.(1/J * g[2] × g[3])
        g²= simplify.(1/J * g[3] × g[1])
        g³= simplify.(1/J * g[1] × g[2])
        g_contra = [g¹, g², g³]
        e_contra = [g/norm(g) for g in g_contra]
        G = simplify.([g[i] ⋅ g[j] for i = 1:3, j = 1:3])
        # invG = simplify.(inv(G))
        invG = simplify.([g_contra[i] ⋅ g_contra[j] for i = 1:3, j = 1:3])
        Γ =  simplify.([sum([invG[i,p]*(∂(G[p,j],q[k]) + ∂(G[p,k],q[j]) - ∂(G[j,k],q[p])) for p=1:3])/2 for i=1:3,j=1:3,k=1:3])
        gdet = simplify(det(G))
        @assert J^2==gdet
        return new(q,g,g_contra,e_cov, e_contra, G,invG,J,gdet,Γ)
    end
end


"""
Check if coordinate system is orthogonal
"""
isorthogonal(C::GenericCoordinates) = isdiag(C.G)



## Covariant and contravariant basis
#types:


struct CovariantVector
    r::Vector3D
    C::GenericCoordinates
end

struct ContravariantVector
    r::Vector3D
    C::GenericCoordinates
end

struct CartesianVector
    r::Vector3D
end
#conversion between contravariant and covariant basis and Cartesian

function CovariantVector(x::ContravariantVector)
    return CovariantVector(x.C.invG*x.r,x.C)
end

function ContravariantVector(x::CovariantVector)
    return ContravariantVector(x.C.G*x.r,x.C)
end

CovariantVector(x::CovariantVector) = x
ContravariantVector(x::ContravariantVector) = x

CartesianVector(x::CovariantVector) = CartesianVector(sum([x.r[i]*x.C.g_cov[i] for i=1:3]))
CartesianVector(x::ContravariantVector) = CartesianVector(sum([x.r[i]*x.C.g_contra[i] for i=1:3]))



#algebra for the two bases

function dot(x::CovariantVector,y::CovariantVector)
    @assert y.C.G == x.C.G
    return CovariantVector(dot(x.C.G*x.r,y.r),x.C)
end

function dot(x::ContravariantVector,y::ContravariantVector)
    @assert y.C.G == x.C.G
    return ContravariantVector(dot(x.C.invG*x.r,y.r),x.C)
end

function cross(x::ContravariantVector,y::ContravariantVector)
    @assert y.C.G == x.C.G
    rout = [sum(√x.C.g*[ϵ([i,j,k])*x.r[j]*y.r[k] for j=1:3,k=1:3]) for i=1:3]
    return CovariantVector(rout,y.C)
end

function cross(x::CovariantVector,y::CovariantVector)
    @assert y.C.G == x.C.G
    rout = [sum(1/√x.C.g*[ϵ([i,j,k])*x.r[j]*y.r[k] for j=1:3,k=1:3]) for i=1:3]
    return ContravariantVector(rout,y.C)
end


#calculus

grad(f::Sym,C::CoordinateSystem) = ContravariantVector([∂(f,C.q[i]) for i=1:3],C)

∇ = grad

div(u::CovariantVector) = 1/√u.C.g*sum([∂(√u.C.g*u.r[i] , u.C.q[i]) for i=1:3])
div(u::ContravariantVector) = 1/√u.C.g*sum([∂(√u.C.g*sum([u.C.invG[i,j]*u.r[i] for i=1:3])[j] , u.C.q[j]) for j=1:3])



end # module
