


function ψp0_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    G = γC(Na,Cox)/sqrt(kb*T)
    vg = Vg/kb/T
    vfb = Vfb/kb/T
    temp = (vg-vfb)*0.5-3*(1+G/sqrt(2))
    return temp + sqrt( temp^2 + 6*(vg-vfb) )
end

function ψp_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    G = γC(Na,Cox)/sqrt(kb*T)
    vg = Vg/kb/T
    vfb = Vfb/kb/T
    δpD = 0.0
    ψ0 = ψp0_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    if vg - vfb > 0
        return 1 -exp(-ψ0) + (sqrt(vg-vfb-1 + exp(-ψ0)+(G/2)^2)-G/2)^2
    else
        return -log(1-ψ0+ ( (vg-vfb-ψ0)/G )^2 )
    end
end

function nq_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    G = γC(Na,Cox)/sqrt(kb*T)
    return 1 + G/2/sqrt(ψp_BSIM6(Vg,Vfb,Na,Cox,T,phif))
end

function q_BSIM6(Vg,Vch,Vfb,Na,Cox,T,phib)
    gam = γC(Na,Cox)/sqrt(kb*T)
    vch = Vch/kb/T
    phib = phib/kb/T/nq_BSIM6(Vg,Vfb,Na,Cox,T,phib)
    psip = ψp_BSIM6(Vg,Vfb,Na,Cox,T,phib)
    T8 = 0.5 * (psip + 1.0 + sqrt((psip - 1.0) * (psip - 1.0) + 0.25 * 2.0 * 2.0))
    sqrtpsip = sqrt(T8)
    T9 = 1.0 + gam / (2.0 * sqrtpsip)
    T0 = (1.0 + (gam / (2.0 * sqrtpsip))) / gam
    T1 = psip - 2.0 * phib - vch
    T2 = T1 - log(4.0 * T0 * sqrtpsip);
    T8 = 0.5 * (T2 - 0.201491 - sqrt(T2 * (T2 + 0.402982) + 2.446562));
    sqrtpsisa = sqrtpsip
    if (T8 <= -68.0) 
        T4 = -100.0
        T5 = 20.0
        if (T8 < T4 - 0.5 * T5)
            T3 = exp(T4)
        else 
            if (T8 > T4 + 0.5 * T5)
                T3 = exp(T8)
            else
                T2 = (T8 - T4) / T5
                T6 = T2 * T2
                T3 = exp(T4 + T5 * ((5.0 / 64.0) + 0.5 * T2 + T6 * ((15.0 / 16.0) - T6 * (1.25 - T6))))
            end
        end
        q = T3 * (1.0 + T1 - T8 - log(2.0 * T0 * (T3 * 2.0 * T0 + 2.0 * sqrtpsisa)))
    else 
        T3 = exp(T8)
        sqrtpsisainv = 1.0 / sqrtpsisa
        T4 = 2.0 * T3 + log(T3 * 2.0 * T0 * (T3 *  2.0 * T0 + 2.0 * sqrtpsisa)) - T1
        T5 = 2.0 + (1.0 / T3) + (T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)
        T3 = T3 - T4 / T5
        T4 = 2.0 * T3 + log(T3 * 2.0 * T0 * (T3 * 2.0 * T0 + 2.0 * sqrtpsisa)) - T1
        T5 = 2.0 + (1.0 / T3) + (T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)
        T6 = ((T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa)) * ((T0 + sqrtpsisainv) / (T0 * T3 + sqrtpsisa))
        T7 = -((1.0 / T3) * (1.0 / T3)) - (1.0 / (sqrtpsisa * sqrtpsisa * sqrtpsisa * (T0 * T3 + sqrtpsisa))) - T6
        q  = T3 - (T4 / T5) * (1.0 + T4 * T7 / (2.0 * T5 * T5))
    end
    return q
end

function id_BSIM6(Vg,Vd,Vs,Vfb,Na,Cox,T,phif)
    qs = q_BSIM6(Vg,Vs,Vfb,Na,Cox,T,phif)
    qd = q_BSIM6(Vg,Vd,Vfb,Na,Cox,T,phif)
    nq = nq_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    return 2*nq*Cox*(kb*T)^2*(qs-qd)*(qs+qd+1)
end
function Id_BSIM6(Vg,Vd,Vs,Vfb,Na,Cox,T,phif)
    qs = q_BSIM6(Vg,Vs,Vfb,Na,Cox,T,phif)
    qd = q_BSIM6(Vg,Vd,Vfb,Na,Cox,T,phif)
    nq = nq_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    return μ_BSIM6(Vg,Vs,Vfb,Na,Cox,T,phif)*2*nq*Cox*(kb*T)^2*(qs-qd)*(qs+qd+1)
end

function μ_BSIM6(Vg,Vs,Vfb,Na,Cox,T,phif)
    U0_i = 67e-3
    UTE_i = -1.5
    UA_i = 0.001
    UC_i = 0
    UD_i = 0.001
    UCS_i = 2.0
    UA1_i = 1e-3
    UC1_i = 0.056e-9
    TRatio = T/300.0
    delTemp = T-300.0
    UD1_i = 0
    UCSTE_i = -4.775e-3
    EU_i = 1.5     
    U0_t   = U0_i * pow(TRatio, UTE_i);
    UA_t   = UA_i * hypersmooth(1.0 + UA1_i * delTemp - 1.0e-6, 1.0e-3);
    UC_t   = UC_i * hypersmooth(1.0 + UC1_i * delTemp - 1.0e-6, 1.0e-3);
    UD_t   = UD_i * pow(TRatio, UD1_i);
    UCS_t  = UCS_i * pow(TRatio, UCSTE_i);
    nq = nq_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    qs = q_BSIM6(Vg,Vs,Vfb,Na,Cox,T,phif)
    phit = kb*T
    psip = ψp_BSIM6(Vg,Vfb,Na,Cox,T,phif)
    n = n_BSIM6(γC(Na,Cox),phif,Vs,phit)
    qis = qis_BSIM6(nq,n,phit,qs)
    qbs = qbs_BSIM6(nq,phit,Vg,Vfb,qs,n,psip)
    Eeffs = Eeff_BSIM6(qbs,qis,1/Cox/3.9/ϵ₀)
    T2 = pow(0.5 * (1.0 + (qis / qbs)), UCS_t);
    T3 = (UA_t + UC_t * 0.0) * pow(Eeffs, EU_i) + UD_t / T2;
    T4 = 1.0 + T3;
    return U0_t/Smooth(T4,1.0,0.0015)
end

mutable struct BSIM6Model <: TransistorModel
    mos::Union{MOSStructure,Nothing}
    W::Number
    L::Number
    μ_0
    VFB
    Nb
    Cox
    T0
    ϕb
end

function BSIM6Model(W,L,mos::MOSStructure)
    BSIM6Model(mos,W,L,67e-3,VFB(mos),abs(C(mos)),Cox(mos),Temperature(mos),ϕb(mos))
end

pow(a,b) = a^b
Eeff_BSIM6(qba,qia,tox,ϵratio=11.7/3.9,η=0.5) = 1e-8*(qba+η*qia)/(ϵratio*tox)
T0_BSIM6(ϕb,Vs,phit) = hypersmooth(2*ϕb+Vs/phit,1e-3)
n_BSIM6(γ,ϕb,Vs,phit) = 1+ γ/(2*sqrt(T0_BSIM6(ϕb,Vs,phit)))
qbs_BSIM6(nq,phit,Vg,Vfb,qs,n,psip) = Smooth(nq*phit * (Vg/phit-Vfb/phit - psip - 2.0 * qs * (n - 1.0)),0,0.1);
qis_BSIM6(nq,n,phit,qs) = 2*n*nq*phit*qs
hypersmooth(x,c)= 0.5*(x+sqrt(x*x+4.0*c*c))
Smooth(x,x0,deltax) =  0.5 * (x + x0 + sqrt((x - x0) * (x - x0) + 0.25 * deltax * deltax))
id(vg,vd,vs,m::BSIM6Model) = id_BSIM6(vg,vd,vs,m.VFB,m.Nb,m.Cox,m.T0,m.ϕb)
ψp(vg,m::BSIM6Model) = ψp_BSIM6(vg,m.VFB,m.Nb,m.Cox,m.T0,m.ϕb)
qi(Vg,V,m::BSIM6Model) = q_BSIM6(Vg,V,m.VFB,m.Nb,m.Cox,m.T0,m.ϕb)
#id(vg,vd,vs,m::BSIM6Model) = id_BSIM6(vg,vd,vs,m.VFB,m.Nb,m.Cox,m.T0,m.ϕb)
Id(vg,vd,vs,m::BSIM6Model) = m.W/m.L*Id_BSIM6(vg,vd,vs,m.VFB,m.Nb,m.Cox,m.T0,m.ϕb)

## Temperature Model


