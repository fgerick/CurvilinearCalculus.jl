using Revise, PyCall, CurvilinearCalculus, LinearAlgebra


z = CurvilinearCalculus.symbols("z", real=true)
a,b,c, H,ϕ = CurvilinearCalculus.@syms a b c H ϕ real=true nonzero=true
q = Vector3D(H,ϕ,z);

assume(H,:real)
assume(-H^2+c^2,:positive)

# global_assumptions=CurvilinearCalculus.global_assumptions
s = √(c^2-H^2)/c

x = a*s*cos(ϕ)
y = b*s*sin(ϕ)


cmap = CoordinateMapping(x,y,z);

CS = GenericCoordinates(cmap,q)

CS.J

isorthogonal(CS)

ψ = CurvilinearCalculus.SymFunction("ψ")

simplify((divergence(1/H*∇(ψ(q...),CS)+1/(3H)*(∇(H,CS)×∇(ψ(H,ϕ),CS))×∇(H,CS))(a=>1,b=>1,c=>1))
