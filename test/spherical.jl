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



differentiate()




# CurvilinearCalculus.simplify.(CS.e_cov[1])
#
# r1 = CovariantVector([1,2,0],CS);
#
# r2 = ContravariantVector(CurvilinearCalculus.Vector3D(1,1,0),CS);
#
# f=CurvilinearCalculus.SymFunction("f")
# g=CurvilinearCalculus.SymFunction("g")
#
# r3 = ContravariantVector([0,0,f(r)],CS)
# r2 = CovariantVector([f(r),0,0],CS)
# CS.g_contra
#
# simplify.(curl(∇(f(r,θ,ϕ),CS)).r)
#
#
# simplify(divergence(curl(r3)))
#
# norm(r1)
# simplify(dot(r3,cross(r1,r3)))
# simplify(dot(CartesianVector(r3),CartesianVector(curl(r3))))
# ContravariantVector(r1)
# simplify(divergence(curl(r3)))
# simplify(divergence(curl(r3)))
#
# simplify(divergence(r2))
