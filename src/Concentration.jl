
function Nc(T,S::SemiconductorMaterial)
    T_0 = S.T0
    S.Nc300*(T/T_0)^(3/2)
end


function Nv(T,S::SemiconductorMaterial)
    T_0 = S.T0
    S.Nv300*(T/T_0)^(3/2)
end

function ni(T,f::DistributionFunc,S::SemiconductorMaterial)
    sqrt(n(T,f,S)*p(T,f,S))
end

function N(T,f::DistributionFunc,S::SemiconductorMaterial) 
    return N(T,Ef(T,f,S),f,S)
end

function N(T,E_f,f::DistributionFunc,S::SemiconductorMaterial)
    E_v = Ev(T,S)
    return n(T,E_f,f,S)+Na(S.Dopping,T,E_f,0,E_v)
end

function Pt(T,f::DistributionFunc,S::SemiconductorMaterial)
    return Pt(T,Ef(T,f,S),f,S)
end
function Pt(T,E_f,f::DistributionFunc,S::SemiconductorMaterial)
    E_c = Ec(T,S)
    return p(T,E_f,f,S)+Na(S.Dopping,T,E_f,E_c,0.0)
end

function Nt(T,f::DistributionFunc,S::SemiconductorMaterial)
    return Nt(T,Ef(T,f,S),f,S)
end
function Nt(T,E_f,f::DistributionFunc,S::SemiconductorMaterial)
    E_c = Ec(T,S)
    return n(T,E_f,f,S)+Nd(S.Dopping,T,E_f,E_c,0)
end

function n(T,f::DistributionFunc,S::SemiconductorMaterial)
    return n(T,Ef(T,f,S),Ec(T,S),f,S)
end

function n(T,E_f,f::DistributionFunc,S::SemiconductorMaterial) 
    return n(T,E_f,Ec(T,S),f,S)
end

function n(T,E_f,E_C,f::DistributionFunc,S::SemiconductorMaterial) 
    ΔE = E_f - E_C
    return Nc(T,S)*f(ΔE,T)
end

function P(T,f::DistributionFunc,S::SemiconductorMaterial)
    η = Ef(T,f,S)
    return P(T,η,f,S)
end

function P(T,E_f,f::DistributionFunc,S::SemiconductorMaterial) 
    E_c = Ec(T,S)
    return p(T,E_f,f,S)+Nd(S.Dopping,T,E_f,E_c,0)
end


function p(T,f::DistributionFunc,S::SemiconductorMaterial)
    η = Ef(T,f,S)
    return p(T,η,f,S)
end

function p(T,E_f,f::DistributionFunc,S::SemiconductorMaterial)  
    p(T,E_f,Ev(T,S),f,S)  
end

function p(T,E_f,E_V,f::DistributionFunc,S::SemiconductorMaterial)  
    ΔE = E_V-E_f
    return Nv(T,S)*f(ΔE,T)
end


