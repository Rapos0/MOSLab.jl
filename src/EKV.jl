function Vth_EKV(Vfb,Na,Cox,T,phif)
    GAMMA = γC(Na,Cox)
    PHI = 2*phif
    VTO = Vfb+PHI+GAMMA*sqrt(PHI)
end

function Vth_EKV(mos::MOSStructure)
    Vth_EKV(VFB(mos),abs(C(mos)),Cox(mos),Temperature(mos),ϕb(mos))
end

function Vp_EKV(Vg,Vfb,Na,Cox,T,phif,sigma =0.0 )
    GAMMA = γC(Na,Cox)
    PHI = 2*phif
    VTO = Vfb+PHI+GAMMA*sqrt(PHI)
    VGprime = Vg - VTO + PHI + GAMMA * sqrt(PHI)
    VP = VGprime - PHI - GAMMA* (sqrt(VGprime+(GAMMA/2.0)*(GAMMA/2.0))-(GAMMA/2.0));
    return VP
end

function Vp_EKV(Vg,mos::MOSStructure )
    Vp_EKV(Vg,VFB(mos),abs(C(mos)),Cox(mos),Temperature(mos),ϕb(mos))
end

function n_EKV(Vg,Vfb,Na,Cox,T,phif,sigma =0.0 )
    GAMMA = γC(Na,Cox)
    PHI = 2*phif
    vt = kb*T
    VP = Vp_EKV(Vg,Vfb,Na,Cox,T,phif)
    n = 1.0 + GAMMA / (2.0*sqrt(PHI + VP + 4.0*vt));
    return n
end

function n_EKV(Vg,mos::MOSStructure )
    n_EKV(Vg,VFB(mos),abs(C(mos)),Cox(mos),Temperature(mos),ϕb(mos))
end

function G_EKV(Na,Cox,T)
    γC(Na,Cox)/sqrt(kb*T)
end


function G_EKV(mos::MOSStructure)
    γC(abs(C(mos)),Cox(mos))/sqrt(kb*Temperature(mos))
end

function q_EKV(Vg,V,Vfb,Na,Cox,T,phif,sigma =0.0 )
    VP = Vp_EKV(Vg,Vfb,Na,Cox,T,phif)
    x=(VP-V)/kb/T 
    iff = (log(1.0+exp( x /2.0)))*(log(1.0+exp( x /2.0)));
    return iff
end

function q_EKV(Vg,V,mos::MOSStructure)
   q_EKV(Vg,V,VFB(mos),abs(C(mos)),Cox(mos),Temperature(mos),ϕb(mos))
end

function id_EKV(Vg,Vd,Vs,Vfb,Na,Cox,T,phif,sigma =0.0 )
    pt = kb*T
    Ispec = 2*n_EKV(Vg,Vfb,Na,Cox,T,phif)*pt*pt*Cox
    return Ispec*(q_EKV(Vg,Vs,Vfb,Na,Cox,T,phif,sigma)-q_EKV(Vg,Vd,Vfb,Na,Cox,T,phif,sigma))

end

function Id_EKV(Vg,Vd,Vs,Vfb,Na,Cox,T,phif,sigma =0.0,Kp=5e-2 )
    pt = kb*T
    Ispec = 2*n_EKV(Vg,Vfb,Na,Cox,T,phif,sigma)*pt*pt*Cox
    Kp = Kp*(T/300)^(-1.5)
    return Kp*Ispec*(q_EKV(Vg,Vs,Vfb,Na,Cox,T,phif,sigma)-q_EKV(Vg,Vd,Vfb,Na,Cox,T,phif,sigma))

end


struct EKV26Model <: TransistorModel
    mos::MOSStructure
    W::Number
    L::Number
    μ_0
end

id(vg,vd,vs,m::EKV26Model) = id_EKV(vg,vd,vs,VFB(m.mos),abs(C(m.mos)),Cox(m.mos),Temperature(m.mos),ϕb(m.mos))
Id(vg,vd,vs,m::EKV26Model,sigma=0.0,Kp=5e-2) = Id_EKV(vg,vd,vs,VFB(m.mos),abs(C(m.mos)),Cox(m.mos),Temperature(m.mos),ϕb(m.mos),sigma,Kp)
Vp(Vg,m::EKV26Model) = Vp_EKV(Vg,m.mos)