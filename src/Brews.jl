struct BrewsModel <: TransistorModel
    mos::MOSStructure
    W::Number
    L::Number
    μ
end

function ψsV(Vgs,V,mdl::BrewsModel)
    ϕt = Temperature(mdl.mos)*kb
    f(ψ) = VFB(mdl.mos)-Vgs+ψ+γ(mdl.mos)*sqrt(ψ/ϕt+(ni(mdl.mos)/C(mdl.mos))^2*exp((ψ-V)/ϕt))
    return fzero(f,V)
end

function Id(Vg,Vd,Vs,Vb,mdl::BrewsModel)
    Vgs = Vg-Vs
    ϕt = Temperature(mdl.mos)*kb
    ψsD = ψsV(Vgs,Vd,mdl)
    ψsS = ψsV(Vgs,Vs,mdl)
    G = γ(mdl.mos)
    f(ψ) = (Vgs-VFB(mdl.mos)-ψ)-G+2*ϕt*((Vgs-VFB(mdl.mos)-ψ)+ϵ₀*ϵᵣ(mdl.mos.semiconductor)*q*C(mdl.mos))/((Vgs-VFB(mdl.mos)-ψ)-G)
    D = hcubature(x->f(x[1]), [ψsS], [ψsD];maxevals=500,reltol=1e-12)[1]
    return mdl.W/mdl.L*Cox(mdl.mos)*mdl.μ*D
end

Id(Vg,Vd,Vs,m::BrewsModel) = Id(Vg,Vd,Vs,Vs,m::BrewsModel)
Id(Vg,Vd,m::BrewsModel) = Id(Vg,Vd,0.0,m::BrewsModel)
gm(Vg,Vd,Vs,Vb,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,Vb,m),Vg)
gm(Vg,Vd,Vs,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,m),Vg)
gm(Vg,Vd,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,0.0,m),Vg)
gmId(Vg,Vd,m::BrewsModel;order=5) = gm(Vg,Vd,m;order=order)/Id(Vd,Vg,m)
gds(Vg,Vd,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(Vg,x,m),Vd)
gdsId(Vg,Vd,m::BrewsModel;order=5) = gds(Vg,Vd,m;order=order)/Id(Vd,Vg,m)

