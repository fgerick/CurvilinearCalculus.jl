module CurvilinearCalculus

using SymPy, PyCall, LinearAlgebra, StaticArrays, Combinatorics

export CoordinateSystem, GenericCoordinates, isorthogonal,
        Vector3D, CoordinateMapping, Metric, CCVector, CartesianVector,
        PhysicalVector

#sympy:
export @syms, simplify, ∂, applyassumptions

#algebra & calculus
export dot,cross
export ∇, grad, divergence, curl, laplacian, ∇², Δ



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

"""
    struct CartesianVector

DOCSTRING
"""
struct CartesianVector
    r::Vector3D
end


"""
    struct GenericCoordinates <: CoordinateSystem

Defines a (non-)orthogonal coordinate system with fields

#Arguments:
- `q`: Coordinates
- `g_cov`: Covariant basis vectors (in Cartesian coordinates)
- `e_cov`: Covariant unit vectors (in Cartesian coordinates)
- `g_contra`: Contravariant basis vectors (in Cartesian coordinates)
- `e_contra`: Contravariant unit vectors (in Cartesian coordinates)
- `F`: transformation matrix new coordinates to Cartesian
- `G`: Metric (defined by dot products of covariant basis vectors)
- `invG`: Inverse of metric
- `J`: Volume element
- `g`: J^2 or det(G)
- `Γ`: Christoffel symbols Γⁱⱼₖ → Γ[i,j,k]
"""
struct GenericCoordinates <: CoordinateSystem
    q::SVector{3,Coordinate}
    g_cov::Vector{CartesianVector}
    e_cov::Vector{CartesianVector}
    g_contra::Vector{CartesianVector}
    e_contra::Vector{CartesianVector}
    F::Metric
    G::Metric
    invG::Metric
    J::Sym
    g::Sym
    Γ::SArray{Tuple{3,3,3},Sym}
end

"""
    GenericCoordinates(r::CoordinateMapping, q::SVector{3, Coordinate})

DOCSTRING
"""
function GenericCoordinates(r::CoordinateMapping, q::SVector{3,Coordinate})
    g_cov = [CartesianVector(∂.(r,q)) for q in q] #eq (5.13)
    J = simplify(g_cov[1] ⋅ (g_cov[2] × g_cov[3]))
    g_contra1 = simplify(1/J * g_cov[2] × g_cov[3])
    g_contra2 = simplify(1/J * g_cov[3] × g_cov[1])
    g_contra3 = simplify(1/J * g_cov[1] × g_cov[2])
    g_contra = [g_contra1, g_contra2, g_contra3]
    G = [g_cov[i] ⋅ g_cov[j] for i = 1:3, j = 1:3]
    F = [∂(r[i],q[j]) for i=1:3,j=1:3] #transformation
    invG = [g_contra[i] ⋅ g_contra[j] for i = 1:3, j = 1:3]

    e_cov = [g_cov[i]/√G[i,i] for i in 1:3]
    e_contra = [g_contra[i]/√invG[i,i] for i in 1:3]

    Γ =  [sum([invG[i,p]*(∂(G[p,j],q[k]) + ∂(G[p,k],q[j]) - ∂(G[j,k],q[p])) for p=1:3])/2 for i=1:3,j=1:3,k=1:3]
    gdet = det(G)

    return GenericCoordinates(q,g_cov,e_cov,g_contra,e_contra,F, G,invG,J,gdet,Γ)
end

"""
    isorthogonal(C::GenericCoordinates) = begin

Check if coordinate system is orthogonal
"""
isorthogonal(C::GenericCoordinates) = isdiag(C.G)



"""
    struct PhysicalVector

DOCSTRING
"""
struct PhysicalVector
    r::Vector3D
    C::GenericCoordinates
end

"""
    struct CCVector

DOCSTRING

#Arguments:
- `cov`: DESCRIPTION
- `contra`: DESCRIPTION
- `C`: DESCRIPTION
"""
struct CCVector
    cov::Vector3D #covariant vector components uᵢ, so that 𝐮 = uᵢ𝐠ⁱ
    contra::Vector3D #contravariant vector components uⁱ, so that 𝐮 = uⁱ𝐠ᵢ
    C::GenericCoordinates
end

#use latex present for vectors in notebooks:
import Base.show

Base.show(io::IO, ::MIME"text/latex", x::CCVector) = print(io, SymPy.sympy.latex(Vector(x.cov), mode="equation*"))
Base.show(io::IO, ::MIME"text/latex", x::PhysicalVector) = print(io, SymPy.sympy.latex(Vector(x.r), mode="equation*"))
Base.show(io::IO, ::MIME"text/latex", x::CartesianVector) = print(io, SymPy.sympy.latex(Vector(x.r), mode="equation*"))
Base.show(io::IO, ::MIME"text/latex", x::Vector3D) = print(io, SymPy.sympy.latex(Vector(x), mode="equation*"))

#conversion between contravariant and covariant basis and Cartesian

#transform covariant vector components to contravariant vector components
"""
    cov2contra(covcomps::Vector3D, C::CoordinateSystem) = begin

DOCSTRING
"""
cov2contra(covcomps::Vector3D,C::CoordinateSystem) = C.invG*covcomps #(eq 2.38)

"""
    contra2cov(contracomps::Vector3D, C::CoordinateSystem) = begin

DOCSTRING
"""
contra2cov(contracomps::Vector3D,C::CoordinateSystem) = C.G*contracomps #(eq 2.40)

"""
    CartesianVector(x::CCVector, usecontracomps::Bool=true)

DOCSTRING
"""
function CartesianVector(x::CCVector,usecontracomps::Bool=true)
    usecontracomps ? CartesianVector(sum([x.contra[i]*x.C.g_cov[i] for i=1:3])) : CartesianVector(sum([x.cov[i]*x.C.g_contra[i] for i=1:3]))
end

#alternative transform:
#CartesianVector(x::ContravariantVector) = CartesianVector(CS.F*x.r)

CartesianVector(x::PhysicalVector) = CartesianVector(sum([x.r[i]*x.C.e_cov[i] for i=1:3]))
CartesianVector(x::CartesianVector) = x

PhysicalVector(x::CCVector) = PhysicalVector(x.contra .* .√diag(x.C.G),x.C)


function CCVector(x::PhysicalVector)
    contra = x.r ./ .√diag(x.C.G)
    cov = contra2cov(contra,x.C)
    CCVector(cov,contra,x.C)
end

# ContravariantVector(x::PhysicalVector) = ContravariantVector(x.r ./ .√diag(x.C.G),x.C)
"""
    CCVector(x::CartesianVector, CS::CoordinateSystem)

DOCSTRING
"""
function CCVector(x::CartesianVector,CS::CoordinateSystem)
    contra = Vector3D(inv(CS.F)*x.r)
    cov = contra2cov(contra,CS)
    CCVector(cov,contra,CS)
end
#
# ContravariantVector(x::CartesianVector,CS::CoordinateSystem) = ContravariantVector(inv(CS.F)*x.r,CS)
# CovariantVector(x::CartesianVector,CS::CoordinateSystem) = CovariantVector(ContravariantVector(x,CS))

#arithmetic functions:
import Base.*,Base.-,Base.+,Base./

*(a::Sym,x::PhysicalVector) = PhysicalVector(x.r*a,x.C)
*(x::PhysicalVector,a::Sym) = a*x

*(a::Sym,x::CCVector) = CCVector(x.cov*a,x.contra*a,x.C)
*(x::CCVector,a::Sym) = a*x


/(x::PhysicalVector,a::Sym) = PhysicalVector(x.r ./a,x.C)
/(x::CCVector,a::Sym) = CCVector(x.cov ./a,x.contra ./a, x.C)
#
# +(x::T,y::T) where T<: CCVector = T(x.r .+ y.r, x.C)
# -(x::T,y::T) where T<: CCVector = T(x.r .- y.r, x.C)


+(x::PhysicalVector,y::PhysicalVector) = PhysicalVector(x.r .+ y.r,x.C)
-(x::PhysicalVector,y::PhysicalVector) = PhysicalVector(x.r .- y.r,x.C)

+(x::CCVector,y::CCVector) = CCVector(x.cov .+ y.cov, x.contra .+ y.contra,x.C)



*(a::Sym,x::CartesianVector) = CartesianVector(x.r*a)
*(x::CartesianVector,a::Sym) = a*x

/(x::CartesianVector,a::Sym) = CartesianVector(x.r ./a)

+(x::CartesianVector,y::CartesianVector) = CartesianVector(x.r .+ y.r)
-(x::CartesianVector,y::CartesianVector) = CartesianVector(x.r .- y.r)


*(a::Number,x::T) where T<: CCVector = T(x.cov*a,x.contra*a,x.C)
*(x::T,a::Number) where T<: CCVector = a*x

/(x::T,a::Number) where T<: CCVector = T(x.r ./a,x.C)

# +(x::T,y::T) where T<: CCVector = T(x.r .+ y.r, x.C)
# -(x::T,y::T) where T<: CCVector = T(x.r .- y.r, x.C)

*(a::Number,x::CartesianVector) = CartesianVector(x.r*a)
*(x::CartesianVector,a::Number) = a*x

/(x::CartesianVector,a::Number) = CartesianVector(x.r ./a)

#algebra for the two bases

function dot(x::CCVector,y::CCVector,usecontracomps::Bool=true)
    @assert y.C.G == x.C.G
    usecontracomps ? dot(x.C.G*x.contra,y.contra) : dot(x.C.invG*x.cov,y.cov)
end

# function dot(x::CovariantVector,y::CovariantVector)
#     @assert y.C.G == x.C.G
#     return dot(x.C.invG*x.r,y.r)
# end

dot(x::CartesianVector,y::CartesianVector) = dot(x.r,y.r)

#between different bases:

dot(x::CartesianVector,y::CCVector, args...) = dot(x,CartesianVector(y, args...))
dot(x::CCVector,y::CartesianVector, args...) = dot(y,x, args...)

function cross(x::CCVector,y::CCVector, usecontracomps::Bool=true)
    @assert y.C.G == x.C.G
    cov = Vector3D([sum(x.C.J*[ϵ([i,j,k])*x.contra[j]*y.contra[k] for j=1:3,k=1:3]) for i=1:3]...)
    contra = Vector3D([sum(1/x.C.J*[ϵ([i,j,k])*x.cov[j]*y.cov[k] for j=1:3,k=1:3]) for i=1:3]...)
    CCVector(cov,contra,x.C)
end

cross(x::CartesianVector,y::CartesianVector) = CartesianVector(cross(x.r,y.r))


#calculus (see Aris 1989; p.169-170)

#Aris eq. (7.55.4) covariant differentiation of a contravariant vector
differentiatecontra(A::CCVector,j::Int,k::Int) = ∂(A.contra[j],A.C.q[k]) + sum([A.C.Γ[j,i,k]*A.contra[i] for i=1:3])

#Aris eq. (7.55.5) covariant differentiation of a covariant vector
differentiatecov(A::CCVector,j::Int,k::Int) = ∂(A.cov[j],A.C.q[k]) - sum([A.C.Γ[i,j,k]*A.cov[i] for i=1:3])


# maybe define scalar in coordinates to avoid second argument?
"""
    grad(f::Sym, C::CoordinateSystem)

DOCSTRING
"""
function grad(f::Sym,C::CoordinateSystem)
    cov = Vector3D([∂(f,C.q[i]) for i=1:3])
    contra = cov2contra(cov,C)
    CCVector(cov,contra,C)
end
∇ = grad

"""
    divergence(u::CCVector, usecontracomps::Bool=true)

DOCSTRING
"""
function divergence(u::CCVector, usecontracomps::Bool=true)
    usecontracomps ? sum([differentiatecontra(u,i,i) for i=1:3]) : sum([u.C.invG[i,j]*differentiatecov(u,i,j) for i=1:3,j=1:3])
end

"""
    laplacian(f::Sym, C::CoordinateSystem) = begin

DOCSTRING
"""
laplacian(f::Sym,C::CoordinateSystem) = 1/√C.g * sum([∂(√C.g*sum([C.invG[i,j]*∂(f,C.q[i]) for i=1:3]),C.q[j]) for j=1:3])
Δ  = laplacian
∇² = Δ

"""
    curl(u::CCVector, usecontracomps::Bool=true)

DOCSTRING
"""
function curl(u::CCVector, usecontracomps::Bool=true)
    contra = usecontracomps ? [1/√u.C.g*sum([ϵt[i,j,k]*sum([u.C.G[k,p]*differentiatecontra(u,p,j) for p=1:3]) for j=1:3,k=1:3 ]) for i=1:3] : [1/√u.C.g*sum([ϵt[i,j,k]*differentiatecov(u,k,j) for j=1:3,k=1:3 ]) for i=1:3]
    contra = Vector3D(contra)
    cov = contra2cov(contra,u.C)
    CCVector(cov,contra,u.C)
end


#sympy conveniences
simplify(x::PhysicalVector) = PhysicalVector(simplify.(x.r),x.C)
simplify(x::CCVector) = CCVector(simplify.(x.cov),simplify.(x.contra),x.C)

simplify(x::CartesianVector) = CartesianVector(simplify.(x.r))

"""
    simplify(C::GenericCoordinates)

DOCSTRING
"""
function simplify(C::GenericCoordinates)
    return GenericCoordinates(C.q,
                simplify.(C.g_cov),
                simplify.(C.e_cov),
                simplify.(C.g_contra),
                simplify.(C.e_contra),
                simplify.(C.F),
                simplify.(C.G),
                simplify.(C.invG),
                simplify.(C.J),
                simplify.(C.g),
                simplify.(C.Γ))
end




end # module
