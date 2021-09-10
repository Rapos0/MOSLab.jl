using MOSLab 
using Revise
using Symbolics
using  SparseArrays
#= 
using  SparseArrays
using LinearAlgebra
using NLsolve
using Zygote =#
net = Netlist()

I1 = CircuitComponent("I0",CurrentSource(1))
R1 = CircuitComponent("Rds",Resistor(1))
G = CircuitComponent("gm",TransConductance(1))
g = 1
d = 2
s = 0
addComponent(net,R1,Dict("p"=>d,"n"=>s))
addComponent(net,I1,Dict("p"=>d,"n"=>s))
addComponent(net,G,Dict("ip"=>g,"im"=>s,"op"=>d,"om"=>s))

ckt = Circuit(net)
NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 

M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1000.0e-4,1.0e-4))
netT = Netlist()
addComponent(netT,M1,Dict("g"=>1,"d"=>d,"s"=>s))
addComponent(netT,CircuitComponent("Vg",VoltageSource(0.3)),Dict("p"=>g,"n"=>s))
addComponent(netT,CircuitComponent("Vd",VoltageSource(1.0)),Dict("p"=>d,"n"=>s))
ckt = Circuit(netT)


M,Im,V = CircuitFunction(ckt)
Mf = eval(build_function(M,V)[1])
If = eval(build_function(Im,V)[1])
v₀ = ones(length(V))
M1 = Matrix(Mf(v₀))
Im1 = Vector(If(v₀))
v₁ = M1\Im1
M1 = Matrix(Mf(v₁))
Im1 = Vector(If(v₁))





NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 

M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1000.0e-4,1.0e-4))
M2 = CircuitComponent("M2",ACMModel(MOSf(300.0),1.0e-4,1.0e-4))
R1 = CircuitComponent("R1",Resistor(10e3))
R2 = CircuitComponent("R2",Resistor(2.0))
V1 = CircuitComponent("V1",VoltageSource(0.6))
VDD = CircuitComponent("V_{DD}",VoltageSource(1.8))
s =  @variables s
net = Netlist()

# addComponent(net,R1,Dict("p"=>1,"n"=>2))
# addComponent(net,R2,Dict("p"=>2,"n"=>0))
# addComponent(net,VDD,Dict("p"=>1,"n"=>0))
addComponent(net,M1,Dict("d"=>2,"g"=>1,"s"=>0))
#addComponent(net,M2,Dict("d"=>0,"g"=>3,"s"=>2))
addComponent(net,R1,Dict("p"=>2,"n"=>3))
addComponent(net,V1,Dict("p"=>1,"n"=>0))
addComponent(net,VDD,Dict("p"=>3,"n"=>0))

ckt = Circuit(net)
M,IM,V = CircuitFunction(ckt)
IVg(t) = 0.9+0.9*sin(2*π*1e3*t)
dt = 1e-5
T = 0
Vo = []
Viv = []
v₀ = ones(length(V))
for it = 1:2e2
    
    Vi = Vg(T)
    V1.component.V = Vi
    M,IM,V = CircuitFunction(ckt)
    Mf = eval(build_function(M,V)[1])
    If = eval(build_function(IM,V)[1])
    Mm = Mf(v₀)
    Im = If(v₀)
    v₁ = Mm\Im
    iter = 0
    error = 100
    while iter < 100 && error > 1e-6
        v₀ = v₁
        Mm = Mf(v₀)
        Im = If(v₀)
        v₁ = Mm\Im
        iter += 1
        error = abs(maximum(v₁-v₀))
    end 
    if T == 0
        Vo = v₁
    else
        Vo = hcat(Vo,v₁)
    end
    push!(Viv,Vi)
    T += dt 
end
plot(Vo[2,:])