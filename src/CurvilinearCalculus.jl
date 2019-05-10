__precompile__(false)

module CurvilinearCalculus

using SymPy, PyCall, LinearAlgebra, StaticArrays, Combinatorics

export CoordinateSystem, GenericCoordinates, @coordinates, isorthogonal,
        Vector3D, CoordinateMapping, Metric, CovariantVector, ContravariantVector,
        CartesianVector, PhysicalVector

#sympy:
export Q, assume, global_assumptions, simplify, refine, ∂

#algebra & calculus
export dot,cross
export ∇, grad, divergence, curl, laplacian



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
import SymPy.simplify

# define sympy dot as real dot product only to avoid assumptions:
import LinearAlgebra.dot
import LinearAlgebra.cross
import LinearAlgebra.norm

dot(x::Sym, y::Sym) = sum(x.*y)
dot(x::Vector3D, y::Vector3D) = sum(x.*y)


#make assumptions work (kind of?)
sp = pyimport("sympy")
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



struct CartesianVector
    r::Vector3D
end


"""
`GenericCoordinates(r::CoordinateMapping, q::SVector{3,Coordinate})`

Defines a (non-)orthogonal coordinate system with fields

`q` Coordinates,
`g_cov` Covariant basis vectors (in cartesian coordinates),
`g_contra` Contravariant basis vectors (in cartesian coordinates),
`e_cov` Covariant unit vectors (in cartesian coordinates),
`G` Metric (defined by dot products of covariant basis vectors),
`invG` Inverse of metric,
`J` Volume element,
`g` J^2 or det(G),
`Γ` Christoffel symbols Γⁱⱼₖ → Γ[i,j,k]

"""
struct GenericCoordinates <: CoordinateSystem
    q::SVector{3,Coordinate}
    g_cov::Vector{CartesianVector}
    g_contra::Vector{CartesianVector}
    e_cov::Vector{CartesianVector}
    F::Metric
    G::Metric
    invG::Metric
    J::Sym
    g::Sym
    Γ::SArray{Tuple{3,3,3},Sym}

    function GenericCoordinates(r::CoordinateMapping, q::SVector{3,Coordinate})
        g = [CartesianVector(∂.(r,q)) for q in q]
        J = simplify(g[1] ⋅ (g[2] × g[3]))
        g¹= simplify(1/J * g[2] × g[3])
        g²= simplify(1/J * g[3] × g[1])
        g³= simplify(1/J * g[1] × g[2])
        g_contra = [g¹, g², g³]
        G = simplify.([g[i] ⋅ g[j] for i = 1:3, j = 1:3])
        F = [∂(r[i],q[j]) for i=1:3,j=1:3] #transformation
        # invG = simplify.(inv(G))
        invG = simplify.([g_contra[i] ⋅ g_contra[j] for i = 1:3, j = 1:3])

        e_cov = [simplify(g[i]/√G[i,i]) for i in 1:3]
        # e_contra = [simplify.(g_contra[i]/√invG[i,i]) for i in 1:3]

        Γ =  simplify.([sum([invG[i,p]*(∂(G[p,j],q[k]) + ∂(G[p,k],q[j]) - ∂(G[j,k],q[p])) for p=1:3])/2 for i=1:3,j=1:3,k=1:3])
        gdet = simplify(det(G))
        # @assert J^2==gdet
        return new(q,g,g_contra,e_cov,F, G,invG,J,gdet,Γ)
    end
end


"""
Check if coordinate system is orthogonal
"""
isorthogonal(C::GenericCoordinates) = isdiag(C.G)



#maybe later:

# struct Scalar
#     val::Sym
#     C::GenericCoordinates
# end

## Covariant and contravariant basis
#types:
abstract type CCVector end


struct CovariantVector <: CCVector
    r::Vector3D
    C::GenericCoordinates
end

struct ContravariantVector <: CCVector
    r::Vector3D
    C::GenericCoordinates
end

struct PhysicalVector <: CCVector
    r::Vector3D
    C::GenericCoordinates
end

import Base.show

show(x::T) where T<: CCVector = show(Vector(x.r))


#conversion between contravariant and covariant basis and Cartesian

CovariantVector(x::ContravariantVector) = CovariantVector(x.C.invG*x.r,x.C)

ContravariantVector(x::CovariantVector) = ContravariantVector(x.C.G*x.r,x.C)


CovariantVector(x::CovariantVector) = x
ContravariantVector(x::ContravariantVector) = x

CartesianVector(x::CovariantVector) = CartesianVector(sum([x.r[i]*x.C.g_cov[i] for i=1:3]))
CartesianVector(x::ContravariantVector) = CartesianVector(sum([x.r[i]*x.C.g_contra[i] for i=1:3]))
CartesianVector(x::PhysicalVector) = CartesianVector(sum([x.r[i]*x.C.e_cov[i] for i=1:3]))
CartesianVector(x::CartesianVector) = x

PhysicalVector(x::CovariantVector) = PhysicalVector(x.r .* .√diag(x.C.G),x.C)
PhysicalVector(x::ContravariantVector) = PhysicalVector(CovariantVector(x))#.r .* .√diag(x.C.invG) # PhysicalVector(CovariantVector(x))

CovariantVector(x::PhysicalVector) = CovariantVector(x.r ./ .√diag(x.C.G),x.C)
ContravariantVector(x::PhysicalVector) = ContravariantVector(x.r ./ .√diag(x.C.invG),x.C)

CovariantVector(x::CartesianVector,CS::CoordinateSystem) = CovariantVector(inv(CS.F)*x.r,CS)
ContravariantVector(x::CartesianVector,CS::CoordinateSystem) = ContravariantVector(CovariantVector(x,CS))

#arithmetic functions:
import Base.*,Base.-,Base.+,Base./

*(a::Sym,x::T) where T<: CCVector = T(x.r*a,x.C)
*(x::T,a::Sym) where T<: CCVector = a*x

/(x::T,a::Sym) where T<: CCVector = T(x.r ./a,x.C)

+(x::T,y::T) where T<: CCVector = T(x.r .+ y.r, x.C)
-(x::T,y::T) where T<: CCVector = T(x.r .- y.r, x.C)

*(a::Sym,x::CartesianVector) = CartesianVector(x.r*a)
*(x::CartesianVector,a::Sym) = a*x

/(x::CartesianVector,a::Sym) = CartesianVector(x.r ./a)

+(x::CartesianVector,y::CartesianVector) = CartesianVector(x.r .+ y.r)
-(x::CartesianVector,y::CartesianVector) = CartesianVector(x.r .- y.r)



#algebra for the two bases

function dot(x::CovariantVector,y::CovariantVector)
    @assert y.C.G == x.C.G
    return dot(x.C.G*x.r,y.r)
end

function dot(x::ContravariantVector,y::ContravariantVector)
    @assert y.C.G == x.C.G
    return dot(x.C.invG*x.r,y.r)
end

dot(x::CartesianVector,y::CartesianVector) = dot(x.r,y.r)

#between different bases:
dot(x::ContravariantVector,y::CovariantVector) = dot(x,ContravariantVector(y))
dot(x::CovariantVector,y::ContravariantVector) = dot(y,x)
dot(x::CartesianVector,y::CovariantVector) = dot(x,CartesianVector(y))
dot(x::CartesianVector,y::ContravariantVector) = dot(x,CartesianVector(y))
dot(x::CovariantVector,y::CartesianVector) = dot(y,x)
dot(x::ContravariantVector,y::CartesianVector) = dot(y,x)




function cross(x::ContravariantVector,y::ContravariantVector)
    @assert y.C.G == x.C.G
    rout = Vector3D([sum(1/√x.C.g*[ϵ([i,j,k])*x.r[j]*y.r[k] for j=1:3,k=1:3]) for i=1:3]...)
    return CovariantVector(rout,y.C)
end

function cross(x::CovariantVector,y::CovariantVector)
    @assert y.C.G == x.C.G
    rout = Vector3D([sum(√x.C.g*[ϵ([i,j,k])*x.r[j]*y.r[k] for j=1:3,k=1:3]) for i=1:3]...)
    return ContravariantVector(rout,y.C)
end

#between bases, choose output basis as the first input
cross(x::CovariantVector,y::ContravariantVector) = cross(x,CovariantVector(y))
cross(x::ContravariantVector,y::CovariantVector) = cross(x,ContravariantVector(y))

cross(x::CartesianVector,y::CartesianVector) = CartesianVector(cross(x.r,y.r))


#calculus (see Aris 1989; p.169-170)

#Aris eq. (7.55.4) covariant differentiation of a contravariant vector -mismatch ?
differentiate(A::CovariantVector,j::Int,k::Int) = ∂(A.r[j],A.C.q[k]) + sum([A.C.Γ[j,i,k]*A.r[i] for i=1:3])

#Aris eq. (7.55.5) covariant differentiation of a covariant vector
differentiate(A::ContravariantVector,j::Int,k::Int) = ∂(A.r[j],A.C.q[k]) - sum([A.C.Γ[i,j,k]*A.r[i] for i=1:3])


# maybe define scalar in coordinates to avoid second argument?
grad(f::Sym,C::CoordinateSystem) = ContravariantVector([∂(f,C.q[i]) for i=1:3],C)
∇ = grad

# divergence(u::ContravariantVector) = sum([differentiate(u,i,i) for i=1:3])
# divergence(u::CovariantVector) = sum([u.C.invG[i,j]*differentiate(u,i,j) for i=1:3,j=1:3])

divergence(u::CovariantVector) = 1/√u.C.g*sum([∂(√u.C.g*u.r[i] , u.C.q[i]) for i=1:3])
divergence(u::ContravariantVector) = divergence(CovariantVector(u))
# divergence(u::ContravariantVector) = 1/√u.C.g*sum([∂(√u.C.g*sum([u.C.invG[i,j]*u.r[i] for i=1:3]) , u.C.q[j]) for j=1:3])

laplacian(f::Sym,C::CoordinateSystem) = 1/√C.g * sum([∂(√C.g*sum([C.invG[i,j]*∂(f,C.q[i]) for i=1:3]),C.q[j]) for j=1:3])
Δ  = laplacian
∇² = Δ


#
curl(u::CovariantVector) = CovariantVector( [1/√u.C.g*sum([ϵt[i,j,k]*sum([u.C.G[k,p]*differentiate(u,p,j) for p=1:3])
                                                    for j=1:3,k=1:3 ]) for i=1:3] ,u.C)
#
#use equation 5.118 from curvilinear pdf lecture notes Brannon - Curvilinear Analysis in a Euclidean Space

# function curl(u::ContravariantVector)
#     b1 = -1/(√u.C.g)*( differentiate(u,2,3)-differentiate(u,3,2) )
#     b2 = -1/(√u.C.g)*( differentiate(u,3,1)-differentiate(u,1,3) )
#     b3 = -1/(√u.C.g)*( differentiate(u,1,2)-differentiate(u,2,1) )
#     return CovariantVector([b1,b2,b3],u.C)
# end
# curl(u::CovariantVector) = curl(ContravariantVector(u))


curl(u::ContravariantVector) = CovariantVector( [1/√u.C.g*sum([ϵt[i,j,k]*differentiate(u,k,j) for j=1:3,k=1:3 ]) for i=1:3] ,u.C)


#sympy conveniences
simplify(x::T) where T<:CCVector = T(simplify.(x.r),x.C)
simplify(x::CartesianVector) = CartesianVector(simplify.(x.r))

end # module
