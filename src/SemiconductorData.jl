
struct SemiconductorData <: AbstractMaterial
    material::SemiconductorMaterial
    Ec::Number
    Ev::Number
    φf::Number
    φfₚ::Number
    φfₑ::Number
    ψ::Number
    n::Number
    p::Number
    Na::Number
    Nd::Number
    T::Number
    dist::DistributionFunc
#=     μₙ::AbstractMobilityModel
    μₚ::AbstractMobilityModel =#
end

include("MobilityModels.jl")

SemiconductorData(S::SemiconductorData,ψ,φfₑ=0.0*Eunit,φfₚ=0.0*Eunit) = SemiconductorData(S.T,S.dist,S.material;psi=ψ,φfₑ=φfₑ,φfₚ=φfₚ)

function SemiconductorData(T,f::DistributionFunc,S::SemiconductorMaterial;psi=0.0*Eunit,φfₑ=0.0*Eunit,φfₚ=0.0*Eunit)
    E_f = Ef(T,f,S)
    _Ec = Ec(T,S)-E_f-psi
    _Ev = Ev(T,S)-E_f-psi
    _φf = Ei(T,f,S)-E_f
    _φfₑ = φfₑ
    _φfₚ = φfₚ
    _n = n(T,E_f,f,S)
    _p = p(T,E_f,f,S)
    _Na = Na(S.Dopping,T,E_f,0,Ev(T,S))
    _Nd = Nd(S.Dopping,T,E_f,Ec(T,S),Ev(T,S))
    return SemiconductorData(S,_Ec,_Ev,_φf,_φfₑ,_φfₚ,psi,_n,_p,_Na,_Nd,T,f)
end

ψ(S::SemiconductorData) = S.ψ
Eg(S::SemiconductorData) = Ec(S)-Ev(S)
Ec(S::SemiconductorData) = S.Ec
Ev(S::SemiconductorData) = S.Ev
Ef(S::SemiconductorData) = 0*Eunit
Ei(S::SemiconductorData) = S.φf

Na(S::SemiconductorData) = S.Na
Nd(S::SemiconductorData) = S.Nd
ni(S::SemiconductorData) = sqrt(S.n*S.p)
n(S::SemiconductorData,ψ=0.0*Eunit,φfₑ=0.0*Eunit) = n(S.T,ψ+φfₑ,Ec(S),S.dist,S.material)
p(S::SemiconductorData,ψ=0.0*Eunit,φfₚ=0.0*Eunit) = p(S.T,ψ+φfₚ,Ev(S),S.dist,S.material)
#n(S::SemiconductorData) = n(S,0.0*Eunit,φfₑ(S))
#p(S::SemiconductorData) = p(S,0.0*Eunit,φfₚ(S))
φfₚ(S::SemiconductorData) = S.φfₚ
φfₑ(S::SemiconductorData) = S.φfₑ
ϵᵣ(S::SemiconductorData) = S.material.ϵₘ
C(S::SemiconductorData) = Nd(S)-Na(S)
ρ(S::SemiconductorData) = q/ϵ₀/ϵᵣ(S) * (p(S) - n(S)+C(S))
ρ(S::SemiconductorData,ψ) = q/ϵ₀/ϵᵣ(S) * (p(S,ψ) - n(S,ψ)+C(S))
DistFun(S::SemiconductorData,E) = S.dist(E)
Temperature(S::SemiconductorData) = S.T
χ(S::SemiconductorData) = S.material.χ
Nv(S::SemiconductorData) = Nv(Temperature(S),S.material)
Nc(S::SemiconductorData) = Nc(Temperature(S),S.material)



function scalling(S::SemiconductorData,L)
    n_i = ni(S)
    ϕₜ = kb*Temperature(S)
    Vscal = ϕₜ 
    Cscal = max(n_i,abs.(C(S)))
    lscal = L
    λ = ϵ₀*11.7*Vscal/(lscal^2*q*Cscal)
    δ = n_i/Cscal
    μscal = 1000
    Ld0 = sqrt(ϵ₀*11.7* ϕₜ/q/abs(Cscal))
    Rscal = ϕₜ*μscal*Cscal/lscal/lscal
    tscal = lscal^2/ϕₜ/μscal
    Jscal = q*Vscal*μscal*Cscal/lscal
    Dscal = μscal*Vscal
    return Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i

end
function Δχ(S::SemiconductorData,δχ)
    M = Δχ(S.material,δχ)
    SN = SemiconductorData(S.T,S.dist,M)
    S.material = M
    S.Eg = SN.Eg
    S.Ec = SN.Ec
    S.Ev = SN.Ev
    S.Ef = SN.Ef
    S.φf = SN.φf
    S.ψ = SN.ψ
    S.Ei = SN.Ei
    S.n = SN.n
    S.p = SN.p
    S.Na = SN.Na
    S.Nd = SN.Nd
    S.N = SN.N
    S.P = SN.P
    S.T = SN.T
end

function ξ(SD::SemiconductorData,psi)
    return Kingston(SD,psi,SD.dist)
end
function ξ(SD::SemiconductorData,psi,phi)
    return Kingston(SD,psi,phi,SD.dist)
end

function KingstonBKP(SD::SemiconductorData,psi::Number,f::BoltzmanDist)
    #convert to psi definition
    psi = psi +(Ef(SD)-Ei(SD))
    n_i = ni(SD)
    E_d = Ed(SD.material.Dopping)
    if Nd(SD.material.Dopping) != 0
        E_d = Ec(SD)-E_d
    else
        E_d = Ev(SD)+E_d
    end
    ϕt = kb*SD.T
    ωDI = (E_d-Ei(SD))/ϕt 
    uc = Ec(SD)/ϕt 
    ub = (Ef(SD)-Ei(SD))/ϕt 
    us = (psi)/ϕt
    p = (ub-us)*sinh(ub)-(cosh(ub)-cosh(us))
    LD = sqrt(ϵ₀*SD.material.ϵₘ*kb*SD.T/2/q/n_i)
    return ϕt*sqrt(2)/LD*sqrt(abs(p))
end

#= function Kingston(SD::SemiconductorData,psi::Number,f::BoltzmanDist)
    #convert to psi definition
    u = psi 
    n_i = ni(SD)
    phi = kb*SD.T
    epsilon = ϵ₀*11.7
    N_a = Na(SD)   
    Ld = sqrt(2*phi*q*N_a/epsilon)
    F(x) = sqrt(abs(exp(-x/phi)+x/phi-1 + n_i^2/N_a^2*(exp(x/phi)-x/phi-1)))
    return  sign(psi)*Ld*F(u)#*sqrt(abs(q * (N_a ^ 2 * exp(-u / phi) * phi + exp(u / phi) * n_i ^ 2 * phi + (u - phi) * N_a ^ 2 - n_i ^ 2 * (u + phi)) / epsilon / N_a))
end =#


function KingstonBKP(SD::SemiconductorData,psi::Number,f::DistributionFunc)
    E_d = Ed(SD.material.Dopping)
    ϕt = kb*SD.T
    us = psi/ϕt
    si = 1.0
    if Nd(SD.material.Dopping) != 0
        E_d = Ec(SD)-E_d
        si=-1
    else
        E_d = Ev(SD)+E_d
    end    
    Fnd(x) = -ϕt * log1p(2 * exp((Ef(SD) - Ec(SD) + x + E_d) / ϕt)) + ϕt * ((Ef(SD) - Ec(SD) + x + E_d) / ϕt)
    Fna(x) = ϕt * log1p(4 * exp((Ev(SD) - x + E_d - Ef(SD)) / ϕt)) - ϕt * ((Ev(SD) - x + E_d - Ef(SD)) / ϕt)
    p1 = Nd(SD.material.Dopping) * (Fnd(-psi)-Fnd(0))
    p2 = -Na(SD.material.Dopping) * (Fna(psi)-Fna(0))
    ηₚ = (Ev(SD)-Ef(SD))/ϕt
    Pp = ϕt*Nv(SD)*(Primitive(f,ηₚ-us)-Primitive(f,ηₚ))
    ηₙ = (Ef(SD)-Ec(SD))/ϕt
    Pn = -ϕt*Nc(SD)*(Primitive(f,ηₙ+us)-Primitive(f,ηₙ))
    K = q/ϵ₀/ϵᵣ(SD) *(p1+p2+Pp+Pn)   
    return si*sign(psi)*sqrt(abs(K))
end 

function Kingston(SD::SemiconductorData,psi0::Number,f::DistributionFunc)
    E_d = Ed(SD.material.Dopping)
    psi =psi0 + (Ef(SD)-Ei(SD))
    ϕt = kb*SD.T
    si = 1.0
    n_i = ni(SD)
    if Nd(SD.material.Dopping) != 0
        E_d = Ec(SD)-E_d
        si=-1
    else
        E_d = Ev(SD)+E_d
    end    
   
    ωDI = (E_d-Ei(SD))/ϕt 
    uc = Ec(SD)/ϕt 
    ub = (Ef(SD)-Ei(SD))/ϕt 
    us = (psi)/ϕt

    Fnd(x) = -ϕt * log1p(2 * exp((Ef(SD) - Ec(SD) + x + E_d) / ϕt)) + ϕt * ((Ef(SD) - Ec(SD) + x + E_d) / ϕt)
    Fna(x) = ϕt * log1p(4 * exp((Ev(SD) - x + E_d - Ef(SD)) / ϕt)) - ϕt * ((Ev(SD) - x + E_d - Ef(SD)) / ϕt)
    p1 = Nd(SD.material.Dopping)/n_i * log( (1 + 0.5*exp(ωDI - us))/(1 + 0.5*exp(ωDI - ub)) )
    p2 = Na(SD.material.Dopping)/n_i * log( (1 + 0.25*exp(us - ωDI))/(1 + 0.25*exp(ub - ωDI)) )
    ωVI =  (Ev(SD)-Ei(SD))/ϕt
    p3 = -1/(sqrt(pi)/2*f(ωVI))*(2/3*Primitive(f,ωVI-ub)-2/3*Primitive(f,ωVI-us))
    
    ωCI =  (Ec(SD)-Ei(SD))/ϕt
    p4 = 1/(sqrt(pi)/2*f(-ωCI))*(2/3*Primitive(f,us-ωCI)-2/3*Primitive(f,ub-ωCI))

    LD = sqrt(ϵ₀*SD.material.ϵₘ*kb*SD.T/2/q/n_i)

    return  sign(psi0)*ϕt/LD*sqrt(abs(p1+p2+p3+p4))

end 

function Kingston(SD::SemiconductorData,psi0,phib,f::DistributionFunc)
    E_d = Ed(SD.material.Dopping)
    E_f = Ef(SD)
    psi = psi0 + (E_f-Ei(SD))
    ϕt = kb*SD.T
    si = 1.0
    n_i = ni(SD)
    if Nd(SD.material.Dopping) != 0
        E_d = Ec(SD)-E_d
        si=-1
        phibn = 0.0
        phibp = phib
    else
        E_d = Ev(SD)+E_d
        phibn = phib
        phibp = 0.0
    end    
   
    ωDI = (E_d-Ei(SD))/ϕt 
    uc = Ec(SD)/ϕt 
    ub = (E_f-Ei(SD))/ϕt
    us = (psi)/ϕt

    Fnd(x) = -ϕt * log1p(2 * exp((Ef(SD) - Ec(SD) + x + E_d) / ϕt)) + ϕt * ((Ef(SD) - Ec(SD) + x + E_d) / ϕt)
    Fna(x) = ϕt * log1p(4 * exp((Ev(SD) - x + E_d - Ef(SD)) / ϕt)) - ϕt * ((Ev(SD) - x + E_d - Ef(SD)) / ϕt)
    p1 = Nd(SD.material.Dopping)/n_i * log( (1 + 0.5*exp(ωDI - us))/(1 + 0.5*exp(ωDI - ub)) )
    p2 = Na(SD.material.Dopping)/n_i * log( (1 + 0.25*exp(us - ωDI))/(1 + 0.25*exp(ub - ωDI)) )
    ωVI =  (Ev(SD)-Ei(SD))/ϕt
    p3 = -1/(sqrt(pi)/2*f(ωVI))*(2/3*Primitive(f,ωVI-ub+phibp/ϕt)-2/3*Primitive(f,ωVI-us+phibp/ϕt))
    
    ωCI =  (Ec(SD)-Ei(SD))/ϕt
    p4 = 1/(sqrt(pi)/2*f(-ωCI))*(2/3*Primitive(f,us-ωCI-phibn/ϕt)-2/3*Primitive(f,ub-ωCI-phibn/ϕt))

    LD = sqrt(ϵ₀*SD.material.ϵₘ*kb*SD.T/2/q/n_i)

    return  sign(psi0)*ϕt/LD*sqrt(abs(p1+p2+p3+p4))

end 




