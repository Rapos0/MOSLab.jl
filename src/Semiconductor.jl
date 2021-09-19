
"""
# Arguments 
    ϵₘ the relative permitivity of the Material 
    Nv300 The effective number of states on the valence band at T0 [cm^-3]
    Nc300 The effective number of states on the conduction band at T0 [cm^-3]
    Eg::BandGapModel Band gap model that allows temperature dependent band gap the BandGap should be in [eV]
    Dopping::Dopping the Ionized dopants models with Dopping on a specific node [cm^-3]
    χ Material Work Function energy [eV]
    T0 Optional parameter to change the temperature in which the Nv and Nc are given [K]
Calculate the total number of acceptors Used in dopping
"""
abstract type AbstractMaterial
end

struct SemiconductorMaterial
    ϵₘ::Number
    Nv300::Number
    Nc300::Number
    Eg::BandGapModel
    Dopping::Dopping
    χ::Number
    T0::Number
end

"""
# Arguments 
    S::SemiconductorMaterial 
Calculate the total number of acceptors Used in dopping in [cm^-3]
"""
Na(S::SemiconductorMaterial) = Na(S.Dopping)
"""
# Arguments 
    S::SemiconductorMaterial 
Calculate the total number of acceptors Used in dopping [cm^-3]
"""
Nd(S::SemiconductorMaterial) = Nd(S.Dopping)

include("Concentration.jl")

"""
# Arguments 
    S::SemiconductorMaterial
    δχ Change on the vacum energy level [eV]
Apply Potential energy difference to the semiconductor vacum energy
"""
function Δχ(S::SemiconductorMaterial,δχ)
    M = SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,S.Dopping,S.χ+δχ,S.T0)
    return M
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    S::SemiconductorMaterial 
Calculate the energy difference between Conduction and the vacum energy
"""
function Ec(T,S::SemiconductorMaterial)
    return -S.χ
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    S::SemiconductorMaterial 
Calculate the energy difference between Valence and the vacum energy
"""
function Ev(T,S::SemiconductorMaterial)
    return -S.χ-Eg(T,S.Eg)
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    S::SemiconductorMaterial 
Intrinsic carriers concentration [cm^-3] using Maxwell-Boltzmann Distribution
"""
function niB(T,S::SemiconductorMaterial)
    return sqrt(Nc(T,S)*Nv(T,S))*exp(-Eg(T,S.Eg)/2.0/kb/T)
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    f:: DistributionFunc The used Distribution Function For carrier concentration calculation
    S::SemiconductorMaterial 
Calculate the Equilibrium Fermi Level [eV] of the semiconductor in relation to the vacum energy
using a Bisection algorithm
"""
function Ef(T,f::DistributionFunc,S::SemiconductorMaterial)
    Ec_S = Ec(T,S)
    Ev_S = Ev(T,S)

    g(x) = N(T,x,f,S)-P(T,x,f,S)
    if sign(g(Ev_S)) != sign(g(Ec_S))
        return find_zero(g,(Ev_S,Ec_S))
    else
        return Na(S.Dopping) > Nd(S.Dopping)  ? Ev_S : Ec_S
    end
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    f:: DistributionFunc The used Distribution Function For carrier concentration calculation
    S::SemiconductorMaterial 
Calculate the Equilibrium Fermi Level [eV] of the semiconductor in relation to the vacum energy
using a Bisection algorithm
"""
function Ef(T,f::BoltzmanDist,S::SemiconductorMaterial)
    E_i = Ei(T,f,S)
    if(Na(S) > 0)
        #P-Type Material
        ΔEf = kb*T*log(Na(S)/niB(T,S))
        E = E_i-ΔEf
    elseif (Nd(S) > 0)
        #N-Type Material
        ΔEf = kb*T*log(Nd(S)/niB(T,S))
        E = E_i+ΔEf
    else
        E = E_i
    end
   return E
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    f:: DistributionFunc The used Distribution Function For carrier concentration calculation
    S::SemiconductorMaterial 
Calculate the Equilibrium Intrinsic Fermi Level [eV] of the semiconductor in relation to the vacum energy
using a Bisection algorithm
"""

function Ei(T,f::DistributionFunc,S::SemiconductorMaterial)
    Ec_s = Ec(T,S)
    Ev_s = Ev(T,S)
    E = find_zero(x-> ustrip( (n(T,x,f,S)-p(T,x,f,S)) ),0.5*(Ev(T,S)+Ec(T,S)))
    if E > Ec_s
        return Ec_s
    elseif E < Ev_s
        return Ev_s
    end
    return E
end

function Ei(T,f::BoltzmanDist,S::SemiconductorMaterial)
    Ec_s = Ec(T,S)
    Ev_s = Ev(T,S)
    E = Ec_s-Eg(T,S.Eg)*0.5-kb*T*log(sqrt(Nc(T,S)/Nv(T,S)))
    return E
end


"""
# Arguments 
    T:: Number Latice Temperature [K]
    f:: DistributionFunc The used Distribution Function For carrier concentration calculation
    S::SemiconductorMaterial 
Calculate the difference between equilibrium Fermi level and The Intrinsic Fermi level [eV]"""
function ϕb(T,f::DistributionFunc,S::SemiconductorMaterial)
    return (Ef(T,f,S)-Ei(T,f,S))/q
end


SetBGModel(S::SemiconductorMaterial,M::BandGapModel) = SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,M,S.Dopping,S.χ,S.T0)

Silicon() = SemiconductorMaterial(11.7,1.04e19/Dunit^3,2.8e19/Dunit^3,VarshniModel(;α=4.9e-4*Eunit/Tunit,β=655*Tunit,Eg0 =1.169*Eunit),Intrinsic(),4.17*Eunit,300*Tunit)
SiliconNU() = SemiconductorMaterial(11.7,1.04e19/Dunit^3,2.8e19/Dunit^3,VarshniModel(;α=4.9e-4*Eunit/Tunit,β=655*Tunit,Eg0 =1.169*Eunit),Intrinsic(),4.17*Eunit,300*Tunit)




