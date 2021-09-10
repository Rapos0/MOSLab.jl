using MOSLab 
using Revise
using Symbolics
using  SparseArrays
using Plots

NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acc.eptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 

M1 = CircuitComponent("M1",ACMModel(MOSf(300.0),1.0,1.0))
M2 = CircuitComponent("M2",ACMModel(MOSf(300.0),1.0,1.0))
#R1 = CircuitComponent("R1",Resistor(100e3))
V1 = CircuitComponent("V1",VoltageSource(0.6))
VDD = CircuitComponent("V_{DD}",VoltageSource(1.8))
net = Netlist()
addComponent(net,M1,Dict("d"=>2,"g"=>1,"s"=>0))
addComponent(net,M2,Dict("d"=>3,"g"=>3,"s"=>2))
addComponent(net,V1,Dict("p"=>1,"n"=>0))
addComponent(net,VDD,Dict("p"=>3,"n"=>0))

ckt = Circuit(net)
M,IM,V = CircuitFunction(ckt)
IVg(t) = 0.9+0.9*sin(2*π*1e3*t)
dt = 1e-5
T = 0
Vo = dc_op(ckt)
#= 
Viv = []
for it = 1:2e2
    Vi = IVg(T)
    V1.component.V = Vi
    if T == 0
        Vo = dc_op(ckt)
    else
        Vo = hcat(Vo,dc_op(ckt))
    end
    push!(Viv,Vi)
    T += dt 
end
plot(Vo[2,:])
plot!(Vo[1,:]) =#