using MOSLab
using Revise


ckt = Netlist()

NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 

M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1.0e-4,1.0e-4))

addComponent(ckt,M1,Dict("d"=>1,"s"=>0,"g"=>2))


cons = Dict("d"=>1,"s"=>0,"g"=>2)
#for c in cons
#    println(last(c) ∈ values(cons))
#end
pos(ckt.nodes[1])