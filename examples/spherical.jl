using Revise, CurvilinearCalculus, LinearAlgebra

# q = CurvilinearCalculus.@syms r θ ϕ

q = CurvilinearCalculus.@coordinates r θ ϕ;
r,θ,ϕ = q;
CurvilinearCalculus.@syms a b c real=true nonzero=true
x = r*cos(ϕ)*sin(θ)
y = r*sin(ϕ)*sin(θ)
z = r*cos(θ)

cmap = CurvilinearCalculus.CoordinateMapping(x,y,z)

CS = CurvilinearCalculus.GenericCoordinates(cmap,q)

CurvilinearCalculus.isorthogonal(CS)




r1 = CurvilinearCalculus.CovariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

r2 = CurvilinearCalculus.ContravariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

CurvilinearCalculus.ContravariantVector(r1)
