using Revise, PyCall, CurvilinearCalculus, LinearAlgebra



# q = CurvilinearCalculus.@coordinates H ϕ z;
# H,ϕ,z = q;
z = CurvilinearCalculus.symbols("z", real=true)
a,b,c, H,ϕ = CurvilinearCalculus.@syms a b c H ϕ real=true nonzero=true
q = Vector3D(H,ϕ,z)
# assume(a,:positive)
# assume(b,:positive)
# assume(c,:positive)
# Q=CurvilinearCalculus.Q
assume(H,:real)
# assume(c,:positive)
# assume(b,:positive)
# assume(a,:positive)

assume(-H^2+c^2,:positive)

# global_assumptions=CurvilinearCalculus.global_assumptions
s = √(c^2-H^2)/c

x = a*s*cos(ϕ)
y = b*s*sin(ϕ)

# CurvilinearCalculus.refine(abs(c^2-H^2))
# global_assumptions

cmap = CurvilinearCalculus.CoordinateMapping(x,y,z)

CS = CurvilinearCalculus.GenericCoordinates(cmap,q)
# CS

CS.J

CurvilinearCalculus.isorthogonal(CS)

r1 = CurvilinearCalculus.CovariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

r2 = CurvilinearCalculus.ContravariantVector(CurvilinearCalculus.Vector3D(1,0,0),CS);

ContravariantVector(r2)
