struct QuadraticModel <: TransistorModel
    mos::MOSStructure
    Vth
    Cox
    W
    L
    μ
end

function Vth_Quadratic(Vfb,Na,tox,T,phif,sigma=0.0)
    phit = kb*T
    G = γ(Na,tox)
    k = 0#log(n/(n-1))
    #uprintln(n,",",phit,",",G,",",k,",",Vfb,",",Na,",",tox,",",T,",",phif)
    return Vfb+2*phif+phit*k+G*sqrt(2*phif)
end

function Vth_Quadratic(mos::MOSStructure)
    return Vth_Quadratic(VFB(mos),abs(C(mos)),mos.tox,Temperature(mos),ϕb(mos))
end

QuadraticModel(W,L,mos::MOSStructure;μ=0.520) = QuadraticModel(mos,Vth_Quadratic(mos),Cox(mos),W,L,μ)


function Id(Vg,Vd,Vs,Vb,mdl::QuadraticModel)
    Vgs = Vg-Vs
    Vds = Vd-Vs
    ϕt = Temperature(mdl.mos)*kb
    @unpack mos,Vth,Cox,W,L,μ = mdl
    Vov = Vgs-Vth
    S = W/L
    if Vov < 0.0
        return S*μ*Cox*ϕt*(1+exp(Vov/ϕt))*(exp(Vds/ϕt))
    elseif Vds < Vov
        return S*μ*Cox*(Vov*Vds-Vds^2*0.5)
    else
        return S*μ*Cox*0.5*(Vov)^2
    end

end

Id(Vg,Vd,Vs,m::QuadraticModel) = Id(Vg,Vd,Vs,Vs,m::QuadraticModel)
# Symbolics.@register Id(Vg,Vd,Vs,m::QuadraticModel)
Id(Vg,Vd,m::QuadraticModel) = Id(Vg,Vd,0.0,m::QuadraticModel)
gm(Vg,Vd,Vs,Vb,m::QuadraticModel) = Zygote.ForwardDiff.derivative(x->Id(x,Vd,Vs,Vb,m),Vg)
IdNR(vg,vd,vs,m::QuadraticModel) = Id(vg,vd,vs,m)-gm(vg,vd,vs,m)*(vg-vs)-gds(vg,vd,vs,m)*(vd-vs)
# Symbolics.@register IdNR(vg,vd,vs,m::QuadraticModel)
# Symbolics.@register gm(Vg,Vd,Vs,Vb,m::QuadraticModel)  
gm(Vg,Vd,Vs,m::QuadraticModel) = Zygote.ForwardDiff.derivative(x->Id(x,Vd,Vs,m),Vg)
gm(Vg,Vd,m::QuadraticModel) = Zygote.ForwardDiff.derivative(x->Id(x,Vd,0.0,m),Vg)
gmId(Vg,Vd,m::QuadraticModel) = gm(Vg,Vd,m)/Id(Vd,Vg,m)
gds(Vg,Vd,m::QuadraticModel) = Zygote.ForwardDiff.derivative(x->Id(Vg,x,m),Vd)
# ModelingToolkit.@register gds(Vg,Vd,Vs,Vb,m::QuadraticModel)  

gdsId(Vg,Vd,m::QuadraticModel;order=5) = gds(Vg,Vd,m)/Id(Vd,Vg,m)
