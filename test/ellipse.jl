@syms z real=true
@syms a b c H s ϕ real=true positive=true
@syms m n integer=true
q = Vector3D(s,ϕ,z);

ψ,ξ,A = CurvilinearCalculus.SymFunction("ψ,ξ,A");

∂=diff
∫=integrate
H = √(1-s^2)
xc=a*s*cos(ϕ)
yc=b*s*sin(ϕ)
zc = c*z
cmap = CoordinateMapping(xc,yc,zc)
CS = GenericCoordinates(cmap,q)

@test !isorthogonal(CS)


G = Metric([a^2*cos(ϕ)^2+b^2*sin(ϕ)^2 (b^2-a^2)*s*sin(ϕ)*cos(ϕ) 0
            (b^2-a^2)*s*sin(ϕ)*cos(ϕ) s^2*(a^2*sin(ϕ)^2+b^2*cos(ϕ)^2) 0
            0 0 c^2])

Ginv = inv(G)

@test simplify.(G)==simplify.(CS.G)
@test simplify.(Ginv)==simplify.(CS.invG)


γ = a*b*c*s
@test CS.J == γ

g_contra1 = CartesianVector([a*cos(ϕ),b*sin(ϕ),0])
g_contra2 = CartesianVector([-a*s*sin(ϕ),b*s*cos(ϕ),0])
g_contra3 = CartesianVector([0,0,c])
g_contra = [g_contra1,g_contra2,g_contra3]


[@test simplify(CS.g_contra[i]).r==simplify(g_contra[i]).r for i=1:3]

F = Metric([a*cos(ϕ) -a*s*sin(ϕ) 0
            b*sin(ϕ) b*s*cos(ϕ) 0
            0 0 c])
invF = Metric([cos(ϕ)/a sin(ϕ)/b 0
                -sin(ϕ)/(a*s) cos(ϕ)/(b*s) 0
                0 0 1/c])

@test CS.F==F
@test simplify.(inv(CS.F))==invF

g_cov1 = CartesianVector([cos(ϕ)/a,sin(ϕ)/b,0])
g_cov2 = CartesianVector([-sin(ϕ)/a,cos(ϕ)/b,0]/s)
g_cov3 = CartesianVector([0,0,1/c])

g_cov = [g_cov1,g_cov2,g_cov3]
Ginv2 = [simplify(dot(g_cov[i],g_cov[j])) for i=1:3,j=1:3]

[@test simplify(CS.g_cov[i]).r==simplify(g_cov[i]).r for i=1:3]

@test simplify.(Ginv)==Ginv2
@test simplify.(Ginv2)==simplify.(CS.invG)


@syms u_x u_y u_z

u = CartesianVector([u_x,u_y,u_z])
ucov=simplify(CovariantVector(u,CS))
ucontra=simplify(ContravariantVector(u,CS))

@test simplify(ucontra).r == simplify.([cos(ϕ)/a*u_x + sin(ϕ)/b*u_y, -sin(ϕ)/(a*s)*u_x + cos(ϕ)/(b*s)*u_y, u_z/c])
# @test simplify(ucontra.r[1]) == cos(ϕ)/a*u_x+sin(ϕ)/b*u_y
# @test simplify(ucontra.r[2]) == -sin(ϕ)/(a*s)*u_x+cos(ϕ)/(b*s)*u_y
# @test simplify(ucontra.r[3]) == u_z/c


@test simplify(CartesianVector(ucov)).r==u.r
@test simplify(CartesianVector(ucontra)).r==u.r

@test simplify(CartesianVector(ucov)).r == simplify(CartesianVector(ucontra)).r

@test simplify(sum([ucontra.r[i]*g_contra[i] for i=1:3])).r==simplify(sum([ucov.r[i]*g_cov[i] for i=1:3])).r



@test ucov.r == simplify(CovariantVector(ucontra)).r
@test ucontra.r == simplify(ContravariantVector(ucov)).r
