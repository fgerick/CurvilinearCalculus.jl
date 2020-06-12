using Revise, CurvilinearCalculus, LinearAlgebra

# q = CurvilinearCalculus.@syms r θ ϕ
@syms r θ ϕ real=true positive=true
q = Vector3D(r,θ,ϕ);


xc = r*cos(ϕ)*sin(θ)
yc = r*sin(ϕ)*sin(θ)
zc = r*cos(θ)
cmap = CoordinateMapping(xc,yc,zc)
CS = GenericCoordinates(cmap,q)
CS = simplify(CS)

isorthogonal(CS)


@symfuns f real=true
u₁,u₂,u₃=CurvilinearCalculus.SymFunction("u₁,u₂,u₃")
# r3 = CovariantVector([u₁(r,θ,ϕ),u₂(r,θ,ϕ),u₃(r,θ,ϕ)],CS)

# CS.g_contra
diag(CS.G)

√CS.G[3,3]

gradu=simplify.([CurvilinearCalculus.differentiate(r3,j,k) for j=1:3,k=1:3])

gradu[1,3]
simplify.(curl(∇(f(r,θ,ϕ),CS)).r)
simplify(divergence(∇(f(r,θ,ϕ),CS))) == simplify(laplacian(f(r,θ,ϕ),CS))
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
