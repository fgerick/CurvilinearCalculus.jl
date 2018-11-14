using Revise, CurvilinearCalculus, LinearAlgebra

# q = CurvilinearCalculus.@syms r θ ϕ

q = CurvilinearCalculus.@coordinates r θ ϕ;
r,θ,ϕ = q;
CurvilinearCalculus.@syms a b c real=true nonzero=true
x = r*cos(ϕ)*sin(θ)
y = r*sin(ϕ)*sin(θ)
z = r*cos(θ)

cmap = CoordinateMapping(x,y,z)

CS = GenericCoordinates(cmap,q)

isorthogonal(CS)




r1 = CovariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

r2 = ContravariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

ContravariantVector(r1)

f=CurvilinearCalculus.SymFunction("f")

CartesianVector(∇(f(r,θ,ϕ),CS))
