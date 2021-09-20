using Test
using MOSLab
SemiUstrip()
S = Silicon()

F = BoltzmanDist()
T = 300.0*Tunit

niB(T,S)
n(T,F,S)
@test isapprox(Eg(T,S.Eg),1.05; atol=0.1)
@test isapprox(log10(ni(T,F,S)),10.0; atol=0.5)
@test n(T,F,S) ≈ ni(T,F,S)
@test p(T,F,S) ≈ ni(T,F,S)
@test Nv(T,S) == 1.04e19
@test Nc(T,S) == 2.8e19
@test p(T,F,S) ≈ n(T,F,S)
@test Ef(T,F,S) ≈ Ei(T,F,S)
@test isapprox(Ei(T,F,S), 0.5*(Ec(T,S)+Ev(T,S)) + kb*T*0.5*(Nv(T,S)/Nc(T,S));atol=0.1)


PD = PDoppingC(1e15)
@test Na(PD) == 1e15
@test Nd(PD) == 0
ND = NDoppingC(1e15)
@test Nd(ND) == 1e15
@test Na(ND) == 0


PS = PSiliconNU(1e15)
@test Na(PS.Dopping) == 1e15
@test P(T,F,PS) ≈ N(T,F,PS) 
@test p(T,F,PS) ≈ Na(PS.Dopping)
@test Ef(T,F,PS) < Ei(T,F,PS)
P(T,F,PS)
n(T,F,PS)
NS = NSiliconNU(1e15)
@test Nd(NS.Dopping) == 1e15
@test P(T,F,NS)≈N(T,F,NS) 
@test n(T,F,NS) ≈ Nd(NS.Dopping)
@test Ef(T,F,NS) > Ei(T,F,NS)


SDI = SemiconductorData(T,F,Silicon();psi=0.0)
@test isapprox(Eg(SDI),1.05; atol=0.1)
@test isapprox(log10(ni(SDI)),10.0; atol=0.5)
@test n(SDI) ≈ ni(SDI)
@test p(SDI) ≈ ni(SDI)
@test p(SDI) ≈ n(SDI)
@test Ef(SDI) ≈ Ei(SDI)
@test isapprox(Ei(SDI), 0.5*(Ec(SDI)+Ev(SDI)) + kb*T*0.5*(Nv(T,S)/Nc(T,S));atol=0.1)