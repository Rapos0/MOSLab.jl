using MOSLab
using Plots
using LaTeXStrings
NpolyDoping = 5.5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5.5e17 # Bulk doping concentration in cm^-3
E_a = 0.0# The ground state energy of the acceptor in eV
t_ox = 3.8e-7 # oxide thickness in cm
t_si = 1e-6
CList2 = reshape( range(colorant"blue", stop=colorant"red",length=7), 1, 7 );

MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 
transistor_model(T) = BSIM6Model(MOSf(T),1.0e-4,1.0e-4)
Idf(Vg,Vd,T) = Id(Vg,Vd,0.0,transistor_model(T)) ## Calculate the Current using the ACM Model of a transistor having the parameters from MOSf(T) and W =1.0 um, L = 1.0 um, similar contructors are available for the other models
plot()
for Vgv in 0.2:0.2:1.8
    plot!(0:0.05:1.8,x->1e6*Idf(Vgv,x,300.0),label="\$V_{GB} = $(Vgv) V\$") # plot Id, Vd characteristics for different VGB
end
PSPG = plot!(xlabel=L"V_{DS} \ [V]",ylabel=L"I_{DS} \ [\mu A]",legend=:topleft)

plot()
for (i,Tv) in enumerate(LinRange(200,500,5))
    plot!(0:0.05:1.0,x->1e6*Idf(x,0.1,Tv),label="\$T = $(Tv) K\$",linecolor=CList2[i]) # plot Id, Vg characteristics for different Temperatures
end
PSPD = plot!(xlabel=L"V_{GS} \ [V]",ylabel=L"I_{DS} \ [\mu A]",legend=:topleft,yaxis=:log10)
plot(PSPG,PSPD) # plot both graphs side by side 
savefig("TransistorCurves.png")

gmIdf(Vg,Vd,T) = gm(Vg,Vd,0.0,transistor_model(T))/Idf(Vg,Vd,T)
plot()
for (i,Tv) in enumerate(LinRange(200,500,5))
    plot!(0:0.01:1.0,x->gmIdf(x,0.1,Tv),label="\$T = $(Tv) K\$",linecolor=CList2[i]) # plot Id, Vg characteristics for different Temperatures
end
PSPD = plot!(xlabel=L"V_{GS} \ [V]",ylabel=L"I_{DS} \ [\mu A]",legend=:topright)
savefig("gmId.png")