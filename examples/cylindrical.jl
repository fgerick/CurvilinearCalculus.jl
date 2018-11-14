using Revise, CurvilinearCalculus, LinearAlgebra

# q = CurvilinearCalculus.@syms r θ ϕ

z = CurvilinearCalculus.symbols("z", real=true)
a,b,c, r,θ = CurvilinearCalculus.@syms a b c r θ real=true nonzero=true
q = Vector3D(r,θ,z)
x = r*cos(θ)
y = r*sin(θ)
z = z

cmap = CoordinateMapping(x,y,z)

CS = GenericCoordinates(cmap,q)

isorthogonal(CS)


CurvilinearCalculus.simplify.(CS.e_cov[1])

r1 = CovariantVector([1,2,0],CS);

r2 = ContravariantVector(CurvilinearCalculus.Vector3D(1,1,0),CS);

f=CurvilinearCalculus.SymFunction("f")
g=CurvilinearCalculus.SymFunction("g")
u₁,u₂,u₃=CurvilinearCalculus.SymFunction("u₁,u₂,u₃")

r3 = CovariantVector([u₁(r,θ,z),u₂(r,θ,z),u₃(r,θ,z)]./.√diag(CS.G),CS)
r4 = ContravariantVector(r3) #[u₁(r,θ,z),u₂(r,θ,z),u₃(r,θ,z)],CS)
CovariantVector(r3) == r3
# r2 = CovariantVector([f(r),0,0],CS)
# CS.g_contra
# CurvilinearCalculus.differentiate(r4,1,2)
# CurvilinearCalculus.differentiate(r3,1,2)

# r2 = CovariantVector([1,0,θ],CS)

simplify(divergence(r4))
simplify(divergence(r3))
# simplify(divergence(r2))

simplify(divergence(∇(f(r,θ,z),CS))) == simplify(CurvilinearCalculus.laplacian(f(r,θ,z),CS))
simplify(divergence(curl(r3)))
# ∇(f(r,θ,z),CS)

CurvilinearCalculus.differentiate(r3,1,3)

gradu=simplify.([CurvilinearCalculus.differentiate(r3,j,k) for j=1:3,k=1:3])

gradu[2,1]
simplify.(curl(∇(f(r,θ,z),CS)).r)

simplify(r2⋅(curl(r2)))

simplify.(curl(r3).r)[1]

norm(r1)
simplify(dot(r3,cross(r1,r3)))
simplify(dot(CartesianVector(r3),CartesianVector(curl(r3))))
ContravariantVector(r1)
simplify(divergence(curl(r3)))
simplify(divergence(curl(r3)))

simplify(divergence(r2))
