abstract type AbstractMobilityModel

end


abstract type AbstractSaturationVelocity
end

struct ConstantVsat <: AbstractSaturationVelocity
    vsat::Number
end

Vsat(v::ConstantVsat,T) = v.vsat

struct SwatzVsat <: AbstractSaturationVelocity
    α::Number
    θ::Number
    T₀::Number
end

Vsat(v::SwatzVsat,T) = v.α/(1+v.θ*exp(T/v.T₀))

struct ConstantMobility <: AbstractMobilityModel
    μ₀::Number
    sat::AbstractSaturationVelocity
    β₀::Number
    βexp::Number
    T₀::Number
end



β(C::ConstantMobility,T) = C.β₀*(T/C.T₀)^C.βexp
Vsat(C::ConstantMobility,T) = Vsat(C.sat,T)
Gsat(C::ConstantMobility,T,ξ) = (1/(1+(C.μ₀*ξ/Vsat(C,T))^β(C,T)))^(1/β(C,T))
μ(C::ConstantMobility,S::SemiconductorData) = C.μ₀*Gsat(C,Temperature(S),ξ(S,ψ(S)))


struct CaugheyThomasMobility <: AbstractMobilityModel
    μ₁::Number
    μ₂::Number
    α::Number
    β::Number
    γ::Number
    δ::Number
    Ncrit::Number
    sat::AbstractSaturationVelocity
    β₀::Number
    βexp::Number
    T₀::Number
end

μ₀(C::CaugheyThomasMobility,T,N) = C.μ₁*(T/C.T₀)^C.α + (C.μ₂*(T/C.T₀)^C.β - C.μ₁*(T/C.T₀)^C.α)/(1+(T/C.T₀)^C.γ*(N/C.Ncrit)^C.δ)
β(C::CaugheyThomasMobility,T) = C.β₀*(T/300.0)^C.βexp
Vsat(C::CaugheyThomasMobility,T) = Vsat(C.sat,T)
Gsat(C::CaugheyThomasMobility,T,N,ξ) = (1/(1+(μ₀(C,T,N)*ξ/Vsat(C,T))^β(C,T)))^(1/β(C,T))
μ(muM::CaugheyThomasMobility,S::SemiconductorData) = μ₀(muM,Temperature(S),abs(C(S)))*Gsat(muM,Temperature(S),abs(C(S)),abs(ξ(S,ψ(S))))
μ(muM::CaugheyThomasMobility,S::SemiconductorData,ψ) = μ₀(muM,Temperature(S),abs(C(S)))*Gsat(muM,Temperature(S),abs(C(S)),abs(ξ(S,ψ)))
μ(muM::CaugheyThomasMobility,S::SemiconductorData,ψ,ϕ) = μ₀(muM,Temperature(S),abs(C(S)))*Gsat(muM,Temperature(S),abs(C(S)),abs(ξ(S,ψ,ϕ)))

μₙCT() = CaugheyThomasMobility(55.24,1429.23,0.0,-2.3,-3.8,0.73,1.072e17,SwatzVsat(2.4e7,0.8,600),2.0,0.01,300.0)
μₚCT() = CaugheyThomasMobility(49.7,479.37,0,-2.2,-3.7,0.7,1.606e17,SwatzVsat(2.4e7,0.8,600),0.6,0.01,300.0)

struct PSPMobility <: AbstractMobilityModel
    FACTUO::Number
    BETN::Number
    Cox::Number
    T₀::Number
    STBET::Number
    THEMU::Number
    STTHEMU::Number
    MUE::Number
    STMUE::Number
    XCOX::Number
    STXCOX::Number
    CS::Number
    STCS::Number
    THCS::Number
    STTHECS::Number
    EFF0::Number
    FETA::Number
    ETAMUSC::Number
    εsi::Number
end

β(PSPm::PSPMobility,T) = PSPm.FACTUO*PSPm.BETN*PSPm.Cox*(T/PSPm.T₀)^PSPm.STBET
θμ(PSPm::PSPMobility,T) = PSPm.THEMU*(T/PSPm.T₀)^PSPm.STTHEMU
μe(PSPm::PSPMobility,T) = PSPm.MUE*(T/PSPm.T₀)^PSPm.STMUE
xcor(PSPm::PSPMobility,T) = PSPm.XCOX*(T/PSPm.T₀)^PSPm.STXCOX
Cs(PSPm::PSPMobility,T) = PSPm.CS*(T/PSPm.T₀)^PSPm.STCS
θCS(PSPm::PSPMobility,T) = PSPm.THCS*(T/PSPm.T₀)^PSPm.STTHECS
ξff0(PSPm::PSPMobility) = PSPm.EFF0*PSPm.Cox/PSPm.εsi
ημ(PSPm::PSPMobility) = PSPm.ETAMUSC*PSPm.FETA
ξeff(PSPm::PSPMobility,ξ) = ξeff0(PSPm)*ξ



#Gmob(PSPm::PSPMobility,ξ,T) = 1+( μe(PSPm,T)*ξeff(PSPm,ξ) )^θμ(PSPm,T) + Cs(PSPm,T) * 

struct CVTMobilityATLAS <: AbstractMobilityModel
    B::Number
    C::Number
    D::Number
    K::Number
    τ::Number
    γ::Number
    μ0::Number 
    μ1::Number 
    μmax::Number 
    CR::Number 
    CS::Number
    α::Number
    β::Number
    PC::Number
    δ::Number
    E::Number
    A::Number
    AL::Number
    η::Number
end

μAC(ξ,C,T,m::CVTMobilityATLAS) = m.B/(ξ^m.E)+(m.C*C^m.τ)/(T*ξ^m.D)
μsr(ξ,n,p,C,m::CVTMobilityATLAS) = m.δ/(ξ^m.K) #m.δ/(ξ^γ(m,n,p,C))
γ(m::CVTMobilityATLAS,n,p,C) = m.A+m.AL*(n+p)/abs(C^m.η)
μb(C,T,m::CVTMobilityATLAS) = m.μ0*exp(-m.PC/C)+(m.μmax*(T/300.0)^(-m.γ)-m.μ0)/(1+(C/m.CR)^m.α)-m.μ1/(1+(m.CS/C)^m.β)
μ(m::CVTMobilityATLAS,ξ,n,p,C,T) = 1.0/(1.0/μAC(ξ,C,T,m)+1.0/μsr(ξ,n,p,C,m)+1.0/μb(C,T,m))
μ(m::CVTMobilityATLAS,S::SemiconductorData,ψ,ϕ_n=0.0,ϕ_p=0.0) = μ(m,abs(ξ(S,ψ)),n(S,ψ-ϕ_n),p(S,ψ-ϕ_p),abs(C(S)),Temperature(S)) 
μcvn() = CVTMobilityATLAS(4.75e7,1.74e5,0.333,2.0,0.125,2.5,52.2,43.4,1417.0,9.68e16,3.43e20,0.68,2.0,0.0,5.82e14,1.0,2.58,5.85e-21,0.0767)


struct CVTMobility <: AbstractMobilityModel
    B
    C
    μ0
    μmax
    μ1
    γ
    CR
    CS
    α
    β
    δ
end

μAC(ξ,C,T,m::CVTMobility) = m.B*T/(ξ)+m.C*(1/ξ^(1/3.0))/(1/T)
μsr(ξ,n,p,C,m::CVTMobility) = m.δ/(ξ^2.0) #m.δ/(ξ^γ(m,n,p,C))
γ(m::CVTMobility,n,p,C) = m.A+m.AL*(n+p)/abs(C^m.η)
μb(C,T,m::CVTMobility) = m.μ0+(m.μmax*(T/300.0)^(-m.γ)-m.μ0)/(1+(C/m.CR)^m.α)-m.μ1/(1+(m.CS/C)^m.β)
μ(m::CVTMobility,ξ,n,p,C,T) = 1.0/(1.0/μAC(ξ,C,T,m)+1.0/μsr(ξ,n,p,C,m)+1.0/μb(C,T,m))
μ(m::CVTMobility,S::SemiconductorData,ψ,ϕ_n=0.0,ϕ_p=0.0) = μ(m,abs(ξ(S,ψ)),n(S,ψ-ϕ_n),p(S,ψ-ϕ_p),abs(C(S)),Temperature(S)) 
μcvn(C) = CVTMobility(4.75e7,1.75e5*C,52.2,1417,23.4,2.5,96.8e16,3.43e20,0.68,2.0,5.82e14)




struct TaschMobility <: AbstractMobilityModel
    RN#2 
    BETAS#2
    MUB#1150
    TMUB#2.5
    D#3.2e-9
    P1#0.09
    B1#1.75
    P2#4.53e-8
    PE#3.17e-7
    B2#-0.25
    Z11#0.0388
    Z22#1.73e-5
    ESR#2.449e7
    BETA#2
    N2#1.1e21
    N1#2e19
    ALPHA#2.0
end
Γ(m::TaschMobility,Ep) = μeff(m)/(1+(μeff(m)*Ep/m.VSAT)^m.BETAS)^(1/m.BETAS)
μph(m::TaschMobility,T,ξ,n) = 1.0/(m.MUBN*(T/300.0)^(-m.TMUBN))+1.0/(Z(m,T,ξ)/m.DN*Y(m,T,n)*(T/300.0)^0.5)
Z(m::TaschMobility,T,ξ) = m.Z11*(T/300.0)/ξ+m.Z22*ξ^(-1/3.0)
Eff(m::TaschMobility,ξ) = ξ+(m.RN-1)*E0(m)/m.RN
Y(m::TaschMobility,T,n)= m.P1*(T/300.0)^B1 + m.p2*n+m.B2*0.0
μsr(m::TaschMobility,ξ) = (m.ESR/ξ)^m.BETA 
μc(m::TaschMobility,C) = m.N2N*(T/300)^1.5/(C*log(1+γB(m,T))-γB(m,T)/(1+γB(m,T)))
γB(m::TaschMobility,T) = m.N1*(T/300.0)^m.ALPHA