using ModelingToolkit, Plots, DifferentialEquations
using MOSLab
using Symbolics
using Revise


NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 
g = 1
d = 2
s = 0
M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1000.0e-4,1.0e-4))
netT = Netlist()
addComponent(netT,M1,Dict("g"=>1,"d"=>d,"s"=>s))
addComponent(netT,CircuitComponent("Vg",VoltageSource(0.3)),Dict("p"=>g,"n"=>s))
addComponent(netT,CircuitComponent("Vd",VoltageSource(1.0)),Dict("p"=>d,"n"=>s))
ckt = Circuit(netT)
dc_op(ckt)