using Revise, CurvilinearCalculus, LinearAlgebra

# q = CurvilinearCalculus.@syms r θ ϕ

z = CurvilinearCalculus.symbols("z", real=true)
a,b,c, r,θ = CurvilinearCalculus.@syms a b c r θ real=true nonzero=true
q = Vector3D(r,θ,z)
x = r*cos(θ)
y = r*sin(θ)
z = z
assume(r,:positive)
assume(r,:nonzero)
cmap = CoordinateMapping(x,y,z)

CS = GenericCoordinates(cmap,q)

isorthogonal(CS)


f=CurvilinearCalculus.SymFunction("f")
g=CurvilinearCalculus.SymFunction("g")
u₁,u₂,u₃=CurvilinearCalculus.SymFunction("u₁,u₂,u₃");

r_physical = PhysicalVector([u₁(r,θ,z),u₂(r,θ,z),u₃(r,θ,z)],CS)
r_cov = CovariantVector(r_physical)
r_contra = CovariantVector(r_physical)


simplify(refine(divergence(r_cov)(sign(r)=>1)))
simplify(refine(divergence(r_contra)(sign(r)=>1)))

CS.G
PhysicalVector(r_contra).r == r_physical.r

simplify.(refine.(map(x->x(sign(r)=>1),PhysicalVector(curl(r_cov)).r)))

divergence(curl(r_contra))

curl(CovariantVector(∇(f(r,θ,z),CS)))


typeof(r4)
typeof(curl(r3))

∇(f(r,θ,z),CS)
CovariantVector(∇(f(r,θ,z),CS))

simplify.(refine.(PhysicalVector(CovariantVector(∇(f(r,θ,z),CS))).r))




simplify.(refine.(PhysicalVector(∇(f(r,θ,z),CS)).r))




typeof(r4)

simplify.(refine.(PhysicalVector(curl(∇(f(r,θ,z),CS)))))

simplify.(refine.(PhysicalVector(curl(CovariantVector(∇(f(r,θ,z),CS))))))

simplify.(refine.(PhysicalVector(CovariantVector(∇(f(r,θ,z),CS)))))

simplify.(refine.(PhysicalVector(∇(f(r,θ,z),CS))))

CS.invG
curl(∇(f(r,θ,z),CS))

simplify(norm(CS.g_cov[1]))

simplify.(curl(CovariantVector(∇(f(r,θ,z),CS).r,CS)).r)

simplify(refine(divergence(curl(r4)))(sign(r)=>1))




simplify.(refine.(curl(CovariantVector(∇(f(r,θ,z),CS))).r))




simplify(divergence(curl(r3)))
# curl()
# simplify.(curl(CovariantVector(∇(f(r,θ,z),CS))).r)
# simplify(CurvilinearCalculus.laplacian(f(r,θ,z),CS))
# simplify(divergence(curl(r3)))
# ∇(f(r,θ,z),CS)

CurvilinearCalculus.differentiate(r3,3,2)

gradu=simplify.([CurvilinearCalculus.differentiate(r3,j,k) for j=1:3,k=1:3])

# gradu[2,1]
simplify.(curl()


simplify.(curl(∇(f(r,θ,z),CS)).r)
simplify(refine(divergence(curl(r3)))(sign(r)=>1))

simplify(r2⋅(curl(r2)))

simplify.(curl(r3).r)[1]

norm(r1)
simplify(dot(r3,cross(r1,r3)))
simplify(dot(CartesianVector(r3),CartesianVector(curl(r3))))
ContravariantVector(r1)
simplify(divergence(curl(r3)))
simplify(divergence(curl(r3)))

simplify(divergence(r2))
