using MOSLab
using Plots
using LaTeXStrings
T = 300.0 # Temperature Kelvin
N_a = 5e16 # Bulk Doping in cm^-3
t_ox = 4e-7 # oxide thickness in cm
t_si = 2.5e-5 # Silicon thickness
Vgv = 1.0 # Gate Voltage
PB = SemiconductorData(T,BoltzmanDist(),PSilicon(N_a))
MOS = MOSStructure(Alluminium(),SiO2(),PB,t_ox,t_si)
sol = PoissonVoronoi(Vgv,MOS,201)
BandDiagram(sol)
plot!(legend=:topright)
savefig("BandDiagram.png")
PBFermi = SemiconductorData(T,FermiDist(),PSilicon(N_a))
MOSFermi = MOSStructure(Alluminium(),SiO2(),PBFermi,t_ox,t_si)
plot(-2.3:0.05:2.3,x-> ψs(x,MOS),label="Boltzmann Distribution")
plot!(-2.3:0.05:2.3,x-> ψs(x,MOSFermi),label="Fermi Distribution",legend=:topleft)
xlabel!("Vg [V]")
ppsi = ylabel!(L"\psi_s \  [V]")
savefig("psi_comp.png")