using MOSLab
using Symbolics
using  SparseArrays
using LinearAlgebra
using NLsolve
using Zygote
net = Netlist()

NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 

M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1.0e-4,1.0e-4))
M2 = CircuitComponent("M2",ACMModel(MOSf(300.0),1.0e-4,1.0e-4))
R1 = CircuitComponent("R1",Resistor(1.0))
R2 = CircuitComponent("R2",Resistor(2.0))
V1 = CircuitComponent("V1",VoltageSource(2.0))
VDD = CircuitComponent("V_{DD}",VoltageSource(5))
s =  @variables s


# addComponent(net,R1,Dict("p"=>1,"n"=>2))
# addComponent(net,R2,Dict("p"=>2,"n"=>0))
# addComponent(net,VDD,Dict("p"=>1,"n"=>0))
addComponent(net,M1,Dict("d"=>2,"g"=>1,"s"=>0))
#addComponent(net,M2,Dict("d"=>0,"g"=>3,"s"=>2))
addComponent(net,R1,Dict("p"=>2,"n"=>3))
addComponent(net,V1,Dict("p"=>1,"n"=>0))
addComponent(net,VDD,Dict("p"=>3,"n"=>0))

ckt = Circuit(net)

f = circuitFuntion(ckt)

R = nlsolve(f,[1.0,1.0,1.0,1.0,1.0],autodiff = :forward)
#props = CircuitProperties(ckt,1,2)
