#= cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Roots
using Plots
using Latexify
using Unitful
using UnitfulAtomic
using UnitfulEquivalences
using UnitfulRecipes
pgfplotsx()

kC =  -14.0*ustrip(auconvert(1.60217662e-19u"C"/8.85418782e-12u"F/m"/4.0/π/11.7))



m = ustrip(auconvert(9.109e-31u"kg"))
V₀ = ustrip(auconvert(-1e6u"eV"))
a = ustrip(auconvert(5.4e-10u"m"))
b0 = ustrip(auconvert(3.9e-15u"m"))
ħ = ustrip(auconvert(6.626070715e-34/2.0/π * u"J*s"))
E = ustrip(auconvert.(collect(LinRange(-200.0,200.0,1000000)).*u"eV")).*(1+0.0im)
q(E) = sqrt(2.0*m*E/(ħ^2))

N = 5
Ve(x,p) = kC * 1/abs(x-p)
r = -a*2.0:b0*50.0:a*2.0
Pot = similar(r).*0.0
for i in 1:N
    Pot += Ve.(r,((i-1)-(N-1)/2)*a)
end
plot(r,Pot,lw=2.0,dpi=300,legend=false)
ylims!(-25,1.0)
xlims!(-1*a,1*a)
rect(x,x0,b,V0) = V0*(((x-x0) > -b/2) & ((x-x0) < b/2))
Ry = .-rect.(r,-a,10000*b0,25.0).-rect.(r,0.0,10000*b0,25.0).-rect.(r,a,10000*b0,25.0)
plot!(r,Ry ,color=:black,legend=false)
k₂(E) = sqrt(2.0*m*(V₀-E)/(ħ^2))

F(E)=cos(q(E)*a)*cosh(k₂(E)*b0) + (k₂(E)^2-q(E)^2)/(2*k₂(E)*q(E))*sin(q(E)*a)*sinh(k₂(E)*b0)
y = F.(E)
plot(real.(y))
gui()
yvalid = copy(y)
yvalid[ab0s.(real.(y)) .> 1.0] .= NaN64

Cf = ustrip(auconvert(1.0u"eV"))
x = real.(acos.(yvalid))./pi
y = real.(E)./Cf

x = vcat(-x,x,-x,x)
y = vcat(y,y,-y,-y)

plot(x,y)
ylims!(-5,5) =#