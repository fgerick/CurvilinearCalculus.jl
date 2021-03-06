
@syms r θ ϕ real=true positive=true
q = Vector3D(r,θ,ϕ)
x = r*cos(ϕ)*sin(θ)
y = r*sin(ϕ)*sin(θ)
z = r*cos(θ)

cmap = CoordinateMapping(x,y,z)
CS = GenericCoordinates(cmap,q)
CS = simplify(CS)
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


u_phys = PhysicalVector([u₁(q...),u₂(q...),u₃(q...)],CS)
u_cc = CCVector(u_phys)
@test PhysicalVector(u_cc).r == u_phys.r
# u_cov = CovariantVector(u_phys)
# u_contra = CovariantVector(u_phys)

# @test u_phys.r == PhysicalVector(u_cov).r
# @test u_phys.r == PhysicalVector(u_contra).r


assumptions=Dict(abs(sin(θ))=>sin(θ),sign(sin(θ))=>1,sign(r)=>1,abs(r)=>r)


@test simplify(divergence(u_cc)(assumptions...)) ==
      simplify(∂(u_phys.r[1]*r^2,r)/r^2 + 1/(r*sin(θ))*∂(u_phys.r[2]*sin(θ),θ) + 1/(r*sin(θ))*∂(u_phys.r[3],ϕ))



curlU_correct = simplify.([(∂(u_phys.r[3]*sin(θ),θ)-∂(u_phys.r[2],ϕ))/(r*sin(θ)),
         1/r*(∂(u_phys.r[1],ϕ)/sin(θ) - ∂(r*u_phys.r[3],r) ),
         1/r*(∂(r*u_phys.r[2],r) - ∂(u_phys.r[1],θ)  )])
# revise()
curlU = map(x->simplify(x(assumptions...)),PhysicalVector(curl(u_cc)).r)

@test curlU == curlU_correct


gradf = map(x->simplify(x(assumptions...)),PhysicalVector(∇(f(q...),CS)).r)
gradf_correct = [∂(f(q...),r), ∂(f(q...),θ)/r, ∂(f(q...),ϕ)/(r*sin(θ))]

@test gradf == gradf_correct

curlofgrad = map(x->simplify(x(assumptions...)),PhysicalVector(curl(∇(f(q...),CS))).r)

@test curlofgrad == [0,0,0]
@test simplify(divergence(curl(u_cc))(assumptions...)) == 0

@test simplify(laplacian(f(q...),CS)(assumptions...)) == simplify(divergence(∇(f(q...),CS))(assumptions...))
