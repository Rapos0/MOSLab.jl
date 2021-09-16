using MOSLab
using Revise
using Plots
using Surrogates
NpolyDoping = 5.5e17 # Gate Npoly doping concentration in cm^-3
N_a = 5.5e17 # Bulk doping concentration in cm^-3
E_a = 0.044 # The ground state energy of the acceptor in eV
t_ox = 4e-7 # oxide thickness in cm
t_si = 1e-6
MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 
M1 = ACMModel(MOSf(300.0),1e-4,1e-4)
SizeRange = [.8e-4,10e-4]
TRange = [273.15+10,273.15+50]
VoltageRange = [-0.1,1.9]
# M1T = SurrogateTransistor(SizeRange,SizeRange,TRange,VoltageRange,M1)
# plot(0:0.1:1.8,x->Id(x,0.3,0.0,M1))
# scatter(0:0.1:5.0,x->Id(0.1,0.3,x,M1T))
# Id(.2,0.3,0.1,M1T)

function fId(x) 
    setT!(x[3],M1)
    1e9*Id(x[1],x[2],0.0,M1)
end
N_samples = 100
sampling_algorithm = LatinHypercubeSample()
lb = [VoltageRange[1],VoltageRange[1],TRange[1]]
ub = [VoltageRange[2],VoltageRange[2],TRange[2]]
sx = sample(N_samples,lb,ub,sampling_algorithm)
sId = fId.(sx)
suId = Kriging(sx,sId,lb,ub)
plot(0:0.1:1.8,x->suId([x,0.1,300.0]))
opt = SRBF
surrogate_optimize(fId,opt(),lb,ub,suId,sampling_algorithm;maxiters=N_samples)
TV = 0:0.01:1.8
Vgf(x) = suId([x,0.1,300.0])
plot!(TV,Vgf.(TV),ribbon=p->std_error_at_point(suId, [p,0.1,300.0]))

plot(TV,x->Id(0.4,x,0.0,M1))
plot!(TV,x->Id(0.3,x,0.0,M1))
plot!(TV,x->Id(0.2,x,0.0,M1))
plot!(TV,x->Id(0.1,x,0.0,M1))

plot(TV,x->Id(x,0.1,0.0,M1),yaxis=:log10)