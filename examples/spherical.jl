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


CurvilinearCalculus.simplify.(CS.e_cov[1])

r1 = CovariantVector([1,2,0],CS);

r2 = ContravariantVector(CurvilinearCalculus.Vector3D(1,1,0),CS);

f=CurvilinearCalculus.SymFunction("f")
g=CurvilinearCalculus.SymFunction("g")
u₁,u₂,u₃=CurvilinearCalculus.SymFunction("u₁,u₂,u₃")
r3 = CovariantVector([u₁(r,θ,ϕ),u₂(r,θ,ϕ),u₃(r,θ,ϕ)],CS)

# CS.g_contra
diag(CS.G)

√CS.G[3,3]

gradu=simplify.([CurvilinearCalculus.differentiate(r3,j,k) for j=1:3,k=1:3])

gradu[1,3]
simplify.(curl(∇(f(r,θ,ϕ),CS)).r)

[differentiate()]
CS.Γ

simplify(divergence(curl(r3)))

norm(r1)
simplify(dot(r3,cross(r1,r3)))
simplify(dot(CartesianVector(r3),CartesianVector(curl(r3))))
ContravariantVector(r1)
simplify(divergence(curl(r3)))
simplify(divergence(curl(r3)))

simplify(divergence(r2))
