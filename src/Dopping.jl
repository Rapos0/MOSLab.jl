abstract type Dopping
end

struct Intrinsic <: Dopping
end
Nd(D::Intrinsic,T,Ef,Ec,Ev) = 0.0/Dunit^3
Na(D::Intrinsic,T,Ef,Ec,Ev) = 0.0/Dunit^3
Nd(D::Intrinsic) = 0.0/Dunit^3
Na(D::Intrinsic) = 0.0/Dunit^3
Ed(D::Dopping) = 0.0*Eunit
struct NDoppingC <: Dopping
    Nd::Number
end
Nd(D::NDoppingC,T,Ef,Ec,Ev) = D.Nd
Na(D::NDoppingC,T,Ef,Ec,Ev) = 0/Dunit^3
Nd(D::NDoppingC) = D.Nd
Na(D::NDoppingC) = 0/Dunit^3
function NSilicon(N)
    S = Silicon()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,NDoppingC(N),S.χ,S.T0)
end
function NSiliconNU(N)
    S = SiliconNU()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,NDoppingC(N),S.χ,S.T0)
end

struct PDoppingC <: Dopping
    Na::Number
end
Nd(D::PDoppingC,T,Ef,Ec,Ev) = 0.0/Dunit^3
Na(D::PDoppingC,T,Ef,Ec,Ev) = D.Na
Nd(D::PDoppingC) = 0.0/Dunit^3
Na(D::PDoppingC) = D.Na
function PSilicon(N)
    S = Silicon()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,PDoppingC(N),S.χ,S.T0)
end
function PSiliconNU(N)
    S = SiliconNU()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,PDoppingC(N),S.χ,S.T0)
end
struct NDopping <: Dopping
    Nd::Number
    Ed::Number
end
Nd(D::NDopping,T,Ef,Ec,Ev) = D.Nd/(1+2*exp((Ef-(Ec-D.Ed))/kb/T)) 
Na(D::NDopping,T,Ef,Ec,Ev) = 0.0/Dunit^3
Nd(D::NDopping) = D.Nd
Na(D::NDopping) = 0.0/Dunit^3
Ed(D::NDopping) = D.Ed
function NSilicon(Nd,Ed) 
    S = Silicon()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,NDopping(Nd,Ed),S.χ,S.T0)
end
function NSiliconNU(Nd,Ed) 
    S = SiliconNU()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,NDopping(Nd,Ed),S.χ,S.T0)
end
struct PDopping <: Dopping
    Na::Number
    Ea::Number
end
Na(D::PDopping,T,Ef,Ec,Ev) = D.Na/(1+4*exp((Ev+D.Ea-Ef)/kb/T))
Nd(D::PDopping,T,Ef,Ec,Ev) = 0.0/Dunit^3
Na(D::PDopping) = D.Na
Nd(D::PDopping) = 0.0/Dunit^3
Ed(D::PDopping) = D.Ea
function PSilicon(Na,Ea) 
    S = Silicon()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,PDopping(Na,Ea),S.χ,S.T0)
end
function PSiliconNU(Na,Ea) 
    S = SiliconNU()
    SemiconductorMaterial(S.ϵₘ,S.Nv300,S.Nc300,S.Eg,PDopping(Na,Ea),S.χ,S.T0)
end
