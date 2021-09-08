function Qfunc_UICM(x)
    u1 = x + 2.0
    if(u1) < -30
        u1 = -30
        wn1 = 1e-10
    elseif u1 <= 0.6
        x1 = exp(u1-1);
	    wn1 = x1 * (1.0 + 4.0 * x1 / 3.0) / (1.0 + x1 * (7 / 3 + 5 * x1 / 6));
    else
        wn1 = (u1-1) - 24 * (((u1-1) + 2) * (u1-1) - 3) / ((7 * (u1-1) + 58) * (u1-1) + 127);
    end
    zn1 = (u1-1) - wn1 - log(wn1);
    temp3 = 1 + wn1;
    y2 = 2 * temp3 * (temp3 + 2 * zn1/3) - zn1;
    en1 = zn1 * y2 / (temp3 * (y2 - zn1));
    wn1 = wn1 * (1 + en1)
    f(q) = abs(q)-1.0+log(abs(q))-x

    try
        wn1 = fzero(f,wn1)
    catch e
        @warn "q didn't converge for x = $x returning default algorithm value"
    end
    return wn1   
end

function ϕsa_UICM(Vg,Vfb,Na,tox,T)
    G = γ(Na,tox)
    phit = kb*T
    #f(x) = (Vg-Vfb-x)^2/G^2-phit*(exp(-x/phit)-1)-x
    #return fzero(f,Vg-Vfb,Roots.Order0(),atol=1e-6, rtol=1e-6)
    alpha=-G*sqrt(phit);    
    phisa=0;
    ug=Vg-Vfb 
    gammasquare = G*G
    phisa = 0.0
    if (ug<alpha)   
        phisa = - 2 * phit * log( ug / alpha )
    end	                        
    if (ug > phit)  
        auxphisat=1+4*(ug-phit)/(gammasquare)
        phisa=ug+(gammasquare/2)*(1-sqrt(auxphisat))          
    end
    if (ug>0)  
        phisa = ug<phisa ? ug : phisa
	  
	else  
        phisa = phisa>ug ? phisa : ug
    end 
    f(x) = (Vg-Vfb-x)^2/G^2-(phit*(exp(-x/phit)-1)+x)
    return fzero(f,phisa,Roots.Order0(),atol=1e-12,rtol=1e-6,maxevals=10000)

end

function n_UICM(Vg,Vfb,Na,tox,T)
    G = γ(Na,tox)
    phit = kb*T
    phisa =ϕsa_UICM(Vg,Vfb,Na,tox,T)
    signphisa = sign(phisa)
    phisa2= abs(phisa);
    safesqrt(x) = x> 0.0 ? sqrt(x) : 0.0
	if (phisa2<=0.00008)	
	    phisa2=-0.00002008;
	    n=1+(G*(1-exp(-(phisa2)/phit)))/(-2*safesqrt((phisa2)+phit*(exp(-(phisa2)/phit))-phit));
	else 
	    n=1+(G*(1-exp(-(phisa)/phit)))/(2*signphisa*safesqrt((phisa)+phit*(exp(-(phisa)/phit))-phit));
    end
    return n
end

function Vp_UICM(Vg,Vfb,Na,tox,T,phif,sigma =0.0 )
    G = γ(Na,tox)
    phit = kb*T
    phisa = ϕsa_UICM(Vg,Vfb,Na,tox,T)
    η = n_UICM(Vg,Vfb,Na,tox,T)
	vp = phisa - 2 * phif - phit*(1+log(η/(η-1)))
    #vp = vp + sigma / n * (vd + vs);
	return vp
end

function Vth_UICM(Vfb,Na,tox,T,phif,sigma=0.0)
    phit = kb*T
    G = γ(Na,tox)
    k = 0#log(n/(n-1))
    #uprintln(n,",",phit,",",G,",",k,",",Vfb,",",Na,",",tox,",",T,",",phif)
    Vt0 = Vfb+2*phif+phit*k+G*sqrt(2*phif+phit*k)
    n = n_UICM(Vt0,Vfb,Na,tox,T)
    k = log(n/(n-1))
    return Vfb+2*phif+phit*k+G*sqrt(2*phif+phit*k)
end

function id_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    Cox = ϵ₀*3.9/tox
    return Cox/2.0*n*phit^2*((qis^2-2qis)-(qid^2-2*qid))
end

function Cgs_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    #n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    Cox = ϵ₀*3.9/tox
    α = (qid+1)/(qis+1)
    return (2/3.0)*Cox*(1+2*α)/(1+α)/(1+α) * qis/(1+qis)
end

function Cgd_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    #n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    Cox = ϵ₀*3.9/tox
    α = (qid+1)/(qis+1)
    return (2/3.0)*Cox*(α*α+2*α)/(1+α)/(1+α) * qid/(1+qid)
end

function Cbs_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    η = n_UICM(Vg,Vfb,Na,tox,T)
    return (η-1)*Cgs_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma)
end


function Cds_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    η = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    Cox = ϵ₀*3.9/tox
    α = (qid+1)/(qis+1)
    return (-4.0/15.0)*η*Cox*(α+3*α+α*α)/(1+α)/(1+α)/(1+α) *qis/(1+qis)
end

function Cbd_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    η = n_UICM(Vg,Vfb,Na,tox,T)
    return (η-1.0)*Cgd_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma)
end
function Csd_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    Cds_UICM(Vg,Vd,Vs,Vfb,Na,tox,T,phif)
end

function Cgb_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    η = n_UICM(Vg,Vfb,Na,tox,T)
    Cox = ϵ₀*3.9/tox
    return (η-1)/η*(Cox-Cgs_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif)-Cgd_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif))
end

function α_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    return  (qid+1)/(qis+1)
end

function ic_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif,sigma=0.0)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    #n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    #Cox = ϵ₀*3.9/tox
    IF = (qis-1)^2-1.0
    IR = (qid-1)^2-1.0

    return IF-IR
end

function Id_UICM(Vg,Vs,Vd,Vfb,Na,tox,T,phif;sigma=0.0,Kp=5e-2,β=-1.5)
    vp = Vp_UICM(Vg,Vfb,Na,tox,T,phif)
    phit = kb*T
    n = n_UICM(Vg,Vfb,Na,tox,T)
    qis = Qfunc_UICM((vp-Vs)/phit)
    qid = Qfunc_UICM((vp-Vd)/phit)
    Cox = ϵ₀*3.9/tox
    return Kp*(T/300)^(β)*Cox/2.0*n*phit^2*((qis^2-2qis)-(qid^2-2*qid))
end

struct ACMModel <: TransistorModel
    Vfb
    Nb
    tox
    T
    ϕB
    sigma::Number
    W::Number
    L::Number
    μ_0
end

ACMModel(m::MOSStructure,W=1.0,L=1.0;sigma=0.0,μ_0=5e-2) = ACMModel(VFB(m),abs(C(m)),m.tox,Temperature(m),ϕb(m),sigma,W,L,μ_0)

id(vg,vd,vs,m::ACMModel) = id_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
ic(vg,vd,vs,m::ACMModel) = ic_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Id(vg,vd,vs,m::ACMModel;β=-1.5) = Id_UICM(vg,vs,vd,m.Vfb,m.Nb,m.tox,m.T,m.ϕB; sigma = m.sigma, β=β)
Vp(Vg,m::ACMModel) = Vp_UICM(Vg,m.Vfb,m.Nb,m.tox,m.T,m.ϕB)
nq(Vg,m::ACMModel) = n_UICM(Vg,m.Vfb,m.Nb,m.tox,m.T)
ϕsa(Vg,m::ACMModel) = ϕsa_UICM(Vg,m.Vfb,m.Nb,m.tox,m.T)
function If(Vg,Vd,m::ACMModel) 
    vp = Vp(Vg,m)
    phit = kb*m.T
    qid = Qfunc_UICM((vp-Vd)/phit)
    return (qid+1)^2-1.0
end
Vth(m::ACMModel) = Vth_UICM(m.Vfb,m.Nb,m.tox,m.T,m.ϕB)

Cgs(vg,vd,vs,m::ACMModel) = Cgs_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Cgd(vg,vd,vs,m::ACMModel) = Cgd_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Cbs(vg,vd,vs,m::ACMModel) = Cbs_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Cgb(vg,vd,vs,m::ACMModel) = Cgb_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Csd(vg,vd,vs,m::ACMModel) = Csd_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Cds(vg,vd,vs,m::ACMModel) = Cds_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
Cbd(vg,vd,vs,m::ACMModel) = Cbd_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
α(vg,vd,vs,m::ACMModel) = α_UICM(vg,vs,vd,m.Vfb,m.ϕB,m.tox,m.T,m.ϕB,m.sigma)
