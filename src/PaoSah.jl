
mutable struct PaoSahModel <: TransistorModel
    mos::MOSStructure
    W::Number
    L::Number
    μ::AbstractMobilityModel
end

function Fpao(ψ,V,mod::PaoSahModel)
   if Na(mod.mos.semiconductor) > Nd(mod.mos.semiconductor)
    return (abs(ψ-V) < 1e-3) ? 0.0 : -q*(n(mod.mos,ψ-V))/(ξ(mod.mos,ψ))

   else
    return (abs(ψ-V) < 1e-3) ? 0.0 : q*(p(mod.mos,ψ))/(ξ(mod.mos,ψ))
   end
end

function FpaoC(ψ,V,mod::PaoSahModel)
    if Na(mod.mos.semiconductor) > Nd(mod.mos.semiconductor)
        #return μ(m,mod.mos.semiconductor,ψ)*ni(mod.mos)^2/abs(C(mod.mos))*exp((ψ-V)/kb/Temperature(mod.mos))/(ξ(mod.mos,ψ)) 
     return (abs(ψ-V) < 1e-3) ? 0.0 : -μ(mod.μ,mod.mos.semiconductor,ψ,V)* q*(n(mod.mos,ψ-V))/(ξ(mod.mos,ψ))
    else
     return (abs(ψ-V) < 1e-3) ? 0.0 : μ(mod.μ,mod.mos.semiconductor,ψ)*q*(p(mod.mos,ψ-V))/(ξ(mod.mos,ψ))
    end
 end

function id(Vg,Vd,Vs,Vb,mod::PaoSahModel)
    return hcubature(x->Fpao(x[1],x[2],mod), [0.0,Vs-Vb], [ψs(Vg,mod.mos),Vd];maxevals=500,reltol=1e-12)[1]
end

function Id(Vg,Vd,Vs,Vb,mod::PaoSahModel)
    return hcubature(x->FpaoC(x[1],x[2],mod), [0.0,Vs-Vb], [ψs(Vg,mod.mos),Vd];maxevals=500,reltol=1e-12)[1]
end

id(Vg,Vd,Vs,m::PaoSahModel) = id(Vg,Vd,Vs,Vs,m::PaoSahModel)
Id(Vg,Vd,Vs,m::PaoSahModel) = Id(Vg,Vd,Vs,Vs,m::PaoSahModel)
id(Vg,Vd,m::PaoSahModel) = id(Vg,Vd,0.0,m::PaoSahModel)
Id(Vg,Vd,m::PaoSahModel) = Id(Vg,Vd,0.0,m::PaoSahModel)
gm(Vg,Vd,Vs,Vb,m::PaoSahModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,Vb,m),Vg)
gm(Vg,Vd,Vs,m::PaoSahModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,m),Vg)
gm(Vg,Vd,m::PaoSahModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,0.0,m),Vg)
gmId(Vg,Vd,m::PaoSahModel;order=5) = gm(Vg,Vd,m;order=order)/Id(Vd,Vg,m)
gds(Vg,Vd,m::PaoSahModel;order=5) = central_fdm(order,1)(x->Id(Vg,x,m),Vd)
gdsId(Vg,Vd,m::PaoSahModel;order=5) = gds(Vg,Vd,m;order=order)/Id(Vd,Vg,m)