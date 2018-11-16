q = @coordinates r θ ϕ
r,θ,ϕ = q
x = r*cos(ϕ)*sin(θ)
y = r*sin(ϕ)*sin(θ)
z = r*cos(θ)

cmap = CoordinateMapping(x,y,z)
CS = GenericCoordinates(cmap,q)

@test isorthogonal(CS)


Γ = zeros(Sym,3,3,3)
Γ[2,1,2] = Γ[2,2,1] = 1/r
Γ[1,2,2] = -r
Γ[1,3,3] = -r*sin(θ)^2
Γ[3,1,3] = Γ[3,3,1] = 1/r
Γ[3,2,3] = Γ[3,3,2] = 1/tan(θ)
Γ[2,3,3] = -sin(θ)*cos(θ)
Γ = simplify.(Γ)
Γ = SArray{Tuple{3,3,3},Sym}(Γ...)

@test Γ==CS.Γ


f=CurvilinearCalculus.SymFunction("f")
u₁,u₂,u₃=CurvilinearCalculus.SymFunction("u₁,u₂,u₃");

r_physical = PhysicalVector([u₁(q...),u₂(q...),u₃(q...)],CS)
r_cov = CovariantVector(r_physical)
r_contra = CovariantVector(r_physical)

# @test PhysicalVector(∇(f(q...),CS)).r == Vector3D(∂(f(q...),r),∂(f(q...),θ)/r,∂(f(q...),z))

@test simplify.(refine.(curl(∇(f(q...),CS)).r)) == Vector3D(0,0,0)
@test divergence(curl(r_cov)) == 0
@test divergence(curl(r_contra)) == 0

@test simplify(laplacian(f(q...),CS)) == simplify(divergence(∇(f(q...),CS)))
