"""
    struct BrewsModel <: TransistorModel
    # Arguments 
    MOS::Union{MOSStructure,Nothing} Optional the mos Structure that was used to define the other parameters
    W:: Transistor width in cm
    L:: Transistor length
    Cox:: Oxide Thickness
    γ:: Bulk Factor 
    ni:: Intrinsic carriers concentration in cm^-3
    Nb:: Bulk doping concentration in cm^-3
    T0:: Temperature in which parameters are taken
    VFB:: The flatband voltage
    ϵᵣ:: Relative permitivity of the substrate
    μ:: Low Field Mobility 

"""

mutable struct BrewsModel <: TransistorModel
    mos::Union{MOSStructure,Nothing}
    W::Number
    L::Number
    Cox
    γ
    ni
    Nb
    T0
    VFB
    ϵᵣ
    μ
end

BrewsModel(W,L,mos::MOSStructure;μ=0.52) = BrewsModel(mos,W,L,Cox(mos),γ(mos),ni(mos),abs(C(mos)),Temperature(mos),VFB(mos),ϵᵣ(mos),μ)

function ψsV(Vgs,V,mdl::BrewsModel)
    @unpack mos,W,L,Cox,G,n_i,Nb,T0,Vfb,er,μ = mdl
    ϕt = kb*T0
    f(ψ) = Vfb-Vgs+ψ+G*sqrt(ψ/ϕt+(n_i/Nb)^2*exp((ψ-V)/ϕt))
    return fzero(f,V)
end

function Id(Vg,Vd,Vs,Vb,mdl::BrewsModel)
    @unpack mos,W,L,Cox,G,n_i,Nb,T0,VFB,er,μ = mdl
    Vgs = Vg-Vs
    ϕt = T0*kb
    ψsD = ψsV(Vgs,Vd,mdl)
    ψsS = ψsV(Vgs,Vs,mdl)
    f(ψ) = (Vgs-Vfb-ψ)-G+2*ϕt*((Vgs-Vfb-ψ)+ϵ₀*er*q*Nb)/((Vgs-Vfb-ψ)-G)
    D = hcubature(x->f(x[1]), [ψsS], [ψsD];maxevals=500,reltol=1e-12)[1]
    return W/L*Cox*μ*D
end

Id(Vg,Vd,Vs,m::BrewsModel) = Id(Vg,Vd,Vs,Vs,m::BrewsModel)
Id(Vg,Vd,m::BrewsModel) = Id(Vg,Vd,0.0,m::BrewsModel)
gm(Vg,Vd,Vs,Vb,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,Vb,m),Vg)
gm(Vg,Vd,Vs,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,Vs,m),Vg)
gm(Vg,Vd,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(x,Vd,0.0,m),Vg)
gmId(Vg,Vd,m::BrewsModel;order=5) = gm(Vg,Vd,m;order=order)/Id(Vd,Vg,m)
gds(Vg,Vd,m::BrewsModel;order=5) = central_fdm(order,1)(x->Id(Vg,x,m),Vd)
gdsId(Vg,Vd,m::BrewsModel;order=5) = gds(Vg,Vd,m;order=order)/Id(Vd,Vg,m)

