struct SemiconductorPoint{T}
    Data::SemiconductorData
    p::T
end

SemiconductorPoint1D(Data::SemiconductorData,p::Real=0.0) = SemiconductorPoint(Data,p)

function SemiconductorPoint1D(Data::SemiconductorData,p::Real,ψ,φfₑ=0.0*Eunit,φfₚ=0.0*Eunit) 
    D = SemiconductorData(Data,ψ,φfₑ,φfₚ)
    SemiconductorPoint1D(D,p)
end

function SemiconductorPoint1D(S1D::SemiconductorPoint,p)
    SemiconductorPoint1D(S1D.Data,p)
end

Eg(S::SemiconductorPoint) = Eg(S.Data)
Ec(S::SemiconductorPoint) = Ec(S.Data)
Ev(S::SemiconductorPoint) = Ev(S.Data)
Ef(S::SemiconductorPoint) = Ef(S.Data)
Ei(S::SemiconductorPoint) = Ei(S.Data)
#= n(S::SemiconductorPoint) = n(S.Data)
p(S::SemiconductorPoint) = p(S.Data) =#
n(S::SemiconductorPoint,ψ=0.0*Eunit,φfₑ=0.0*Eunit) = n(S.Data,ψ,φfₑ)
p(S::SemiconductorPoint,ψ=0.0*Eunit,φfₚ=0.0*Eunit) = p(S.Data,ψ,φfₚ)
Na(S::SemiconductorPoint) = Na(S.Data)
Nd(S::SemiconductorPoint) = Nd(S.Data)
N(S::SemiconductorPoint) = N(S.Data)
N(S::SemiconductorPoint,ψ) = N(S.Data,ψ)
P(S::SemiconductorPoint) = P(S.Data)
P(S::SemiconductorPoint,ψ) = P(S.Data,ψ)
Pt(S::SemiconductorPoint) = Pt(S.Data)
Nt(S::SemiconductorPoint) = Nt(S.Data)
ni(S::SemiconductorPoint) = ni(S.Data)
φf(S::SemiconductorPoint) = φf(S.Data)
r(S::SemiconductorPoint) = S.p
ψ(S::SemiconductorPoint) = S.Data.ψ
φfₑ(S::SemiconductorPoint) = S.Data.φfₑ
φfₚ(S::SemiconductorPoint) = S.Data.φfₚ
ϵᵣ(S::SemiconductorPoint) = ϵᵣ(S.Data)
ξ(S::SemiconductorPoint) = ξ(S.Data,ψ(S))
ξ(S::SemiconductorPoint,ψ) = ξ(S.Data,ψ)
ξ(S::SemiconductorPoint,ψ,ϕ) = ξ(S.Data,ψ,ϕ)
ρ(S::SemiconductorPoint) = ρ(S.Data)
ρ(S::SemiconductorPoint,ψ) = ρ(S.Data,ψ)
Δχ(S::SemiconductorPoint,δχ) = Δχ(S.Data,δχ)
χ(S::SemiconductorPoint) = χ(S.Data)
Temperature(S::SemiconductorPoint) = Temperature(S.Data)
Nv(S::SemiconductorPoint) = Nv(S.Data)
Nc(S::SemiconductorPoint) = Nc(S.Data)
DistFun(S::SemiconductorPoint,E) = DistFun(S.Data,E)
C(S::SemiconductorPoint) = C(S.Data)

function applyPsi(S::SemiconductorPoint,ψ) 
    D = S.Data
    _p = S.p
    _φfₚ = φfₚ(S)
    _φfₑ = φfₑ(S)
    return SemiconductorPoint1D(D,_p,ψ,_φfₑ,_φfₚ)
end

function applyEfe(S::SemiconductorPoint,φfₑ) 
    D = S.Data
    _p = S.p
    psi = ψ(S)
    _φfₚ = φfₚ(S)
    return SemiconductorPoint1D(D,_p,psi,φfₑ,_φfₚ)
end

function applyEfp(S::SemiconductorPoint,φfₚ) 
    D = S.Data
    _p = S.p
    psi = ψ(S)
    _φfₑ = φfₑ(S)
    return SemiconductorPoint1D(D,_p,psi,_φfₑ,φfₚ)
end