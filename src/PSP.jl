const se = 4.6051701859880916e+02
const se05 = 2.3025850929940458e+02
const ke = 1.0e-200
const ke05 = 1.0e-100
const keinv = 1.0e200
const ke05inv = 1.0e100

## TODO: Create PSP Struct

PSP_VDSx(VDS) = VDS^2/(sqrt(VDS^2+0.01)+0.1)
PSP_VSBx(VSBs,VDS) = VSBs+0.5*(VDS-PSP_VDSx(VDS))

MINA(x,y,a) = (x+y-sqrt((x-y)^2.0+a))/2.0
χ(y::Number) = y^2/(2.0+y^2)
χp(y::Number) = 4*y/(2+y^2)^2
χpp(y::Number) = (8-12*y^2)/(2+y^2)^3.0
ν(a,c) = a+c
μ1(a,c,τ) = (a+c)^2/τ + c^2/2 - a
σ1(a,c,τ,η) = a*ν(a,c)/(μ1(a,c,τ)+(c^2 /3-a)+c*(a+c)/μ1(a,c,τ)) + η
μ2(a,b,c,τ) = (a+c)^2/τ + c^2/2 - a*b
σ2(a,b,c,τ,η) = a*ν(a,c)/(μ2(a,b,c,τ)+(c^2/3-a*b)+c*(a+c)/μ2(a,b,c,τ)) + η

P3(u) =  (1.0 + (u) * (1.0 + 0.5 * ((u) * (1.0 + (u)/3.0))))
ϕx(φb) = 2*0.95*φb
aϕ(φb) = 2.5e-3*(2*φb)^2
bϕ(φb) = 2.5e-3*(2*φb)^2
ϕxs(φb) = MINA(ϕx(φb)-0.5*sqrt(bϕ(φb)),0,aϕ(φb))
ϕV(φb) = MINA(Vs,Vs+Vd,bϕ(φb))+ϕx(φb)
Vsbs(Vs,φb) = Vs-MINA(ϕV(φb),0,aϕ(φb))+ϕxs(φb)

function expl_high(x)
    if ((x) < se05) 
            res       = exp(x)
    else
            res       =  ke05inv * P3((x) - se05)
    end
end

function expl_low(x)
    if ((x) > -se05) 
            res       =  exp(x)
    else 
            res       = ke05 / P3(-se05 - (x))
    end
    return res
end


function sigma2(a,b,c,tau,eta) 
    nu           =  (a) + (c)
    mutau        =  (nu) * (nu) + (tau) * (0.5 * ((c) * (c)) - (a) * (b))
    y            =  (eta) + (a) * nu * (tau) / (mutau + (nu / mutau) * (tau) * (tau) * (c) * ((c) * (c) * 1/3.0 - (a) * (b)))
    return y
end

function sigma(a,c,tau,eta)
    nu           =  (a) + (c)
    mutau        =  nu * nu + (tau) * (0.5 * ((c) * (c)) - (a))
    y            =  (eta) + (a) * nu * (tau) / (mutau + (nu / mutau) * (tau) * (tau) * (c) * ((c) * (c) * 1/3.0 - (a)))
    return y
end


function γ(N,tox)
    Cox = ϵ₀*3.9/tox
    return sqrt(2*q*11.7*ϵ₀*N)/Cox
end

function γC(N,Cox)
    return sqrt(2*q*11.7*ϵ₀*N)/Cox
end

function SP_psi(SD::SemiconductorData,Vg,Vfb,tox,V=0)
    phit = kb*Temperature(SD)
    xg = (Vg-(Vfb))/phit;
    xn = (2.0*abs(Ei(SD))+V)/phit;
    G = γ(abs(C(SD)),tox)/sqrt(phit);
    xsub = xg + G^2.0/2 -G*(xg-1+G^2.0/4)^(1/2)   
    eta = s(xsub,xn+3,5)
    a = (xg-eta)^2-G.^2*eta+G^2
    c = 2.0*(xg-eta)+G^2
    t = xn-eta+log(a/G^2)
    x0 = eta+sigmaSP(a,c,t)
    D0 = exp(x0-xn)
    p = 2*(xg-x0)+G.^2*(1+D0)
    q = (xg-x0)^2-G^2*(x0+D0-1)
    x = x0+2*q./(p+sqrt(p^2-2*(2-G^2*D0)*q))
    return phit*x
end



function PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vsb=0.0)
    phit = kb*T
    φb *= 2.0
    φb = abs(φb)
    phib = φb
    xg = (Vg-Vfb)/phit
    Gf = γ(Na,tox)/sqrt(phit)
    inv_phit1 = 1/phit
    oneSixth = 1/6.0
    invSqrt2 = 1/sqrt(2.0)
    xi           =  1.0 + Gf * invSqrt2
    inv_xi       =  1.0 /  xi
    Ux           =  Vsb * inv_phit1
    xn         =  phib * inv_phit1 + Ux
    if (xn < se)
        delta     =  exp(-xn)
    else 
        delta     = ke / P3(xn - se)
    end
    margin       =  1.0e-5 * xi
    
    inv_xi = 1/xi
    Gf2 = Gf^2
    inv_Gf2 = 1/Gf2
    if (abs(xg) <= margin)  
        SP_S_temp1 =  inv_xi * inv_xi * oneSixth * invSqrt2
        sp         =  xg * inv_xi * (1.0 + xg * (1.0 - (delta)) * Gf * SP_S_temp1)
    else  
        if (xg < -margin)  
            SP_S_yg     = -xg
            SP_S_ysub   = 1.25 * (SP_S_yg * inv_xi)
            SP_S_eta    = 0.5 * (SP_S_ysub + 10 - sqrt((SP_S_ysub - 6.0) * (SP_S_ysub - 6.0) + 64.0))
            SP_S_temp   = SP_S_yg - SP_S_eta
            SP_S_a      = SP_S_temp * SP_S_temp + Gf2*(SP_S_eta + 1.0)
            SP_S_c      = 2.0 * SP_S_temp - Gf2
            SP_S_tau    = -SP_S_eta + log(SP_S_a * inv_Gf2) 
            SP_S_y0 = sigma(SP_S_a, SP_S_c, SP_S_tau, SP_S_eta) 
            SP_S_delta0 = expl_high(SP_S_y0) 
            SP_S_delta1 = 1.0 / SP_S_delta0
            SP_S_temp   = 1.0 / (2.0 + SP_S_y0 * SP_S_y0)
            SP_S_xi0    = SP_S_y0 * SP_S_y0 * SP_S_temp
            SP_S_xi1    = 4.0 * (SP_S_y0 * SP_S_temp * SP_S_temp)
            SP_S_xi2    = (8.0 * SP_S_temp - 12.0 * SP_S_xi0) * SP_S_temp * SP_S_temp
            SP_S_temp   = SP_S_yg - SP_S_y0
            SP_S_temp1  = (delta) * SP_S_delta1
            SP_S_pC     = 2.0 * SP_S_temp + Gf2 * (SP_S_delta0 - 1.0 - SP_S_temp1 + (delta) * (1.0 - SP_S_xi1))
            SP_S_qC     = SP_S_temp * SP_S_temp - Gf2 * (SP_S_delta0 - SP_S_y0 - 1.0 + SP_S_temp1 + (delta) * (SP_S_y0 - 1.0 - SP_S_xi0))
            SP_S_temp   = 2.0 - Gf2 * (SP_S_delta0 + SP_S_temp1 - (delta) * SP_S_xi2)
            SP_S_temp   = SP_S_pC * SP_S_pC - 2.0 * (SP_S_qC * SP_S_temp)
            sp          = -SP_S_y0 - 2.0 * (SP_S_qC / (SP_S_pC + sqrt(SP_S_temp)))
        else 
            SP_xg1    = 1.0 / (1.25 + Gf * 7.324648775608221e-001)
            SP_S_A_fac= (xi * 1.25 * SP_xg1 - 1.0) * SP_xg1
            SP_S_xbar = xg * inv_xi * (1.0 + SP_S_A_fac * xg)
            SP_S_temp = expl_low(-SP_S_xbar) 
            SP_S_w    = 1.0 - SP_S_temp
            SP_S_x1   = xg + Gf2 * 0.5 - Gf * sqrt(xg + Gf2 * 0.25 - SP_S_w)
            SP_S_bx   = (xn) + 3.0
            SP_S_eta  = MINA(SP_S_x1, SP_S_bx, 5.0) - 0.5 * (SP_S_bx - sqrt(SP_S_bx * SP_S_bx + 5.0))
            SP_S_temp = xg - SP_S_eta
            SP_S_temp1= exp(-SP_S_eta)
            SP_S_temp2= 1.0 / (2.0 + SP_S_eta * SP_S_eta)
            SP_S_xi0  = SP_S_eta * SP_S_eta * SP_S_temp2
            SP_S_xi1  = 4.0 * (SP_S_eta * SP_S_temp2 * SP_S_temp2)
            SP_S_xi2  = (8.0 * SP_S_temp2 - 12.0 * SP_S_xi0) * SP_S_temp2 * SP_S_temp2
            SP_S_a    = max(1.0e-40, SP_S_temp * SP_S_temp - Gf2 * (SP_S_temp1 + SP_S_eta - 1.0 - (delta) * (SP_S_eta + 1.0 + SP_S_xi0)))
            SP_S_b    = 1.0 - 0.5 * (Gf2 * (SP_S_temp1 - (delta) * SP_S_xi2))
            SP_S_c    = 2.0 * SP_S_temp + Gf2 * (1.0 - SP_S_temp1 - (delta) * (1.0 + SP_S_xi1))
            SP_S_tau  = (xn) - SP_S_eta + log(SP_S_a / Gf2)
            SP_S_x0 = sigma2(SP_S_a, SP_S_b, SP_S_c, SP_S_tau, SP_S_eta) 
            if (SP_S_x0 < se05) 
                SP_S_delta0 = exp(SP_S_x0)
                SP_S_delta1 = 1.0 / SP_S_delta0
                SP_S_delta0 = (delta) * SP_S_delta0
            else 
                if (SP_S_x0 > (xn) - se05) 
                    SP_S_delta0 = exp(SP_S_x0 - (xn))
                    SP_S_delta1 = (delta) / SP_S_delta0
                else 
                    SP_S_delta0 = ke05 / P3((xn) - SP_S_x0 - se05)
                    SP_S_delta1 = ke05 / P3(SP_S_x0 - se05)
                end 
            end 
            SP_S_temp   = 1.0 / (2.0 + SP_S_x0 * SP_S_x0)
            SP_S_xi0    = SP_S_x0 * SP_S_x0 * SP_S_temp
            SP_S_xi1    = 4.0 * (SP_S_x0 * SP_S_temp * SP_S_temp)
            SP_S_xi2    = (8.0 * SP_S_temp - 12.0 * SP_S_xi0) * SP_S_temp * SP_S_temp
            SP_S_temp   = xg - SP_S_x0
            SP_S_pC     = 2.0 * SP_S_temp + Gf2 * (1.0 - SP_S_delta1 + SP_S_delta0 - (delta) * (1.0 + SP_S_xi1))
            SP_S_qC     = SP_S_temp * SP_S_temp - Gf2 * (SP_S_delta1 + SP_S_x0 - 1.0 + SP_S_delta0 - (delta) * (SP_S_x0 + 1.0 + SP_S_xi0))
            SP_S_temp   = 2.0 - Gf2 * (SP_S_delta1 + SP_S_delta0 - (delta) * SP_S_xi2)
            SP_S_temp   = SP_S_pC * SP_S_pC - 2.0 * (SP_S_qC * SP_S_temp)
            sp          = SP_S_x0 + 2.0 * (SP_S_qC / (SP_S_pC + sqrt(SP_S_temp)))
        end 
    end
    xs = sp
    Es = 0.0
    Ds = 0.0
    Ps = 0.0
    if xg > 0.0
        Es = exp(-xs)
        Ds = abs(1/Es-xs-1-χ(xs))*delta
        Ps = abs(xs-1+Es)
    end
    xgs = xg-xs <= 0 ? xg-xs : Gf*sqrt(Ds+Ps)
    psiss = phit*xs
    return psiss,Es,Ds,Ps,xgs
end

function PSP_psid(Vg,Na,tox,Vfb,T,φb;Vsb=0.0,Vdse=0.0)
    phit = kb*T
    x1 = 1.25
    φb *= 2.0
    xg = (Vg-Vfb)/phit
    G = γ(Na,tox)/sqrt(phit)
    ξ = 1+G/sqrt(2.0)
    xns = (φb+Vsb)/phit
    xnd = (φb+Vsb+Vdse)/phit
    kds = exp(-Vdse/phit)
    Δns = exp(-xns)
    Δnd = Δns*kds
    xmrg = 1e-5*ξ
    if xg > xmrg
        bx= xnd+3.0
        xg1 = x1+G*sqrt(exp(-x1)+x1-1)
        xbar = xg/ξ*(1+xg*(ξ*x1-xg1)/xg1^2)
        x0 = xg+G^2/2.0-G*sqrt(xg+G^2/4-1+exp(-xbar))   
        η = MINA(x0,bx,5) - (bx-sqrt(bx^2+5.0))/2.0
        a = (xg-η)^2-G^2*(exp(-η)+η-1-Δnd*(η+1+χ(η)))
        b = 1-G^2/2.0*(exp(-η)-Δnd*χpp(η))
        c = 2*(xg-η)+G^2*(1-exp(-η)-Δnd*(1+χp(η)))
        τ = xnd-η+log(abs(a/G^2))
        y0 = σ2(a,b,c,τ,η)
        Δ0 = exp(y0)
        p = 2*(xg-y0)+G^2*(1-1/Δ0+Δnd*(Δ0-1-χp(y0)))
        q = (xg-y0)^2-G^2*(y0+1/Δ0-1+Δnd*(Δ0-y0-1-χ(y0)))
        xd = y0+2*q/(p+sqrt(p^2-2*q*(2-G^2*(1/Δ0+Δnd*(Δ0-χpp(y0))))))
    elseif xg <= xmrg
        xd = xg/ξ*(1+G*xg*(1-Δnd)/(ξ^2*6*sqrt(2)))
    else
        xg1 = x1+G*sqrt(exp(-x1)+x1-1)
        xbar = xg/ξ*(1+xg*(ξ*x1-xg1)/xg1^2)
        x0 = xg+G^2/2.0-G*sqrt(xg+G^2/4-1+exp(-xbar))
        bx = xns+3
        η=MINA(x0,bx,5.0)-(bx-sqrt(bx^2+5))/2.0
        a = (xg-η)^2-G^2*(exp(-η)+η-1-Δns*(η+1+χ(η)))
        b = 1-G^2/2.0*(exp(-η)-Δns*χpp(η))
        c = 2*(xg-η)+G^2*(1-exp(-η)-Δns*(1+χp(η)))
        τ = xns-η+log(a/G^2)
        y0 = σ2(a,b,c,τ,η)
        Δ0 = exp(y0)
        p = 2*(xg-y0)+G^2*(1-1/Δ0+Δns*(Δ0-1-χp(y0)))
        q = (xg-y0)^2-G^2*(y0+1/Δ0-1+Δns*(Δ0-y0-1-χ(y0)))
        xd = y0+2*q/(p+sqrt(p^2-2*q*(2-G^2*(1/Δ0+Δns*(Δ0-χpp(y0))))))
    end
    ψs,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vsb)
    xs = ψs*phit
    xds = xd-xs
    if  xds < 1e-10
        p = 2*xgs+G^2*(1-Es+Δnd*(1/Es-1-χp(xs)))
        q = G^2*(1-kds)*Ds
        ξ = 1-G^2/2.0*(Es+Δnd*(1/Es-χpp(xs)))
        xds = 2*q/(p+sqrt(p^2-4*ξ*q))
        xd = xs+xds
    end
    Ed = exp(-xd)
    Dd = (1/Ed-xd-1-χ(xs))*Δnd
    psisd = phit*xd
    return psisd,Ed,Dd,xds
end


function PSP_VaricapPsi(Vg,Na,tox,Vfb,T,φb)
    #psisd,Ed,Dd,Δψ = PSP_psid(Vg,Na,tox,Vfb,T,φb)
    psiss,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb) 
    return psiss
end


function PSP_ψm(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    psiss,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vs)
    psisd,Ed,Dd,Pd,xds = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vd)
    phit = kb*T
    φb *= 2.0
    G = γ(Na,tox)/sqrt(phit)
    x_s = psiss/phit
    x_d = psisd/phit
    xg = (Vg-Vfb)/phit
    xds = x_d-x_s
    if xg > 0.0
        xm = (x_s+x_d)*0.5
        Em = sqrt(Es*Ed)
        Dbar = (Ds+Dd)*0.5
        Dm = Dbar + xds^2/8.0*(Em-2/G^2)
        Pm = xm-1.0+Em
        xgm = G*sqrt(Dm+Pm)
        
        return xm*phit,Em,Dm,Pm,xds,xgm
    else 
        return 0,0,0,0,0,0
    end
end

PSP_qim(G,phit,Dm,Pm,xgm) = G^2*phit*Dm/(xgm+G/sqrt(Pm))
PSP_α(G,Em,Pm) = 1.0 + G*(1-Em)/2/sqrt(Pm)
PSP_qis(G,Ds,xgs,Ps,phit)= G^2*phit*Ds/(xgs+G*sqrt(Ps))
PSP_αs(G,Es,Ps) = 1 + G*(1-Es)/2/sqrt(Ps)
PSP_qbs(phit,G,Ps) = phit*G*sqrt(Ps)
PSP_μx(Vbsx,T,Xcor=0,stxcor=0) = 1+Xcor*Vdsx/(1+0.2*Xcor*Vbsx)
PSP_Eeff(qbs,qis,cox,esi,η=1/2)=1e-8*cox/esi*(qbs+η*qis)
#PSP_ρs(θr,qis) = θr
PSP_θμ(T,STTHEMU=1.5,THEMU=1.5)= THEMU*(T/300.0)^STTHEMU
PSP_μe(T,MUE=0.5,STMUE=0) = MUE*(T/300.0)^STMUE
PSP_Gmob(μe,Eeff,θμ,Cs,qbs,qis,θcs=2.0)=(1+(μe*Eeff)^θμ+Cs*(qbs/(qis+qbs))^θcs)
PSP_μ(μ0,T,Gmob,U0ST) = μ0*(T/300.0)^U0ST/Gmob
PSP_wsat(qis) = 100*qis/(100+qis)
PSP_θsat(T,Gmob,wsat,THESATG=0,THESAT=1,STTHETASAT=1.0) = THESAT*(T/300.0)^STTHETASAT/Gmob*(1+THESATG*wsat)
PSP_ϕ∞(qis,αs,phit)=qis/αs+phit
PSP_ysat(θsat,ϕ∞) = θsat*ϕ∞/sqrt(2)
PSP_za(ysat)=2/(1+sqrt(1+4*ysat))
PSP_ϕ0(ϕ∞,za,ysat) =ϕ∞*za*(1+0.86*za*ysat*(1-za^2*ysat)/(1+4*za^3*ysat^2))
PSP_asat(xgs,G) = xgs+G^2/2.0
PSP_ϕ2(phit,G,Ds,asat) = phit*0.98*G^2*Ds/(asat+sqrt(asat^2-0.98*G^2*Ds))
PSP_ϕsat(ϕ0,ϕ2) = 2*ϕ0*ϕ2/(ϕ0+ϕ2+sqrt((ϕ0+ϕ2)^2-3.95*ϕ0*ϕ2))

function PSP_θsat(Vg,Vs,Na,tox,Vfb,T,φb,Cox)
    psiss,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vs)
    G = γ(Na,tox)/sqrt(phit)
    phit = kb*T
    qis =  PSP_qis(G,Ds,xgs,Ps,phit)
    αs = PSP_αs(G,Es,Ps)
    qbs = PSP_qbs(phit,G,Ps)
    Eeff = PSP_Eeff(qbs,qis,Cox,ϵ₀*11.7)
    μe = PSP_μe(T)
    θμ = PSP_θμ(T)
    CS = 0
    Gmob = PSP_Gmob(μe,Eff,θμ,CS,qbs,qis)
    wsat = PSP_wsat(qis)
    θsat = PSP_θsat(T,Gmob,wsat)
end

function PSP_Vdsat(Vg,Vs,Na,tox,Vfb,T,φb,Cox)
    psiss,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vs)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    
    qis =  PSP_qis(G,Ds,xgs,Ps,phit)
    αs = PSP_αs(G,Es,Ps)
    qbs = PSP_qbs(phit,G,Ps)
    Eeff = PSP_Eeff(qbs,qis,Cox,ϵ₀*11.7)
    μe = PSP_μe(T)
    θμ = PSP_θμ(T)
    CS = 0
    Gmob = PSP_Gmob(μe,Eeff,θμ,CS,qbs,qis)
    wsat = PSP_wsat(qis)
    θsat = PSP_θsat(T,Gmob,wsat)
    ϕ∞ = PSP_ϕ∞(qis,PSP_αs(G,Es,Ps),phit)
    ysat = PSP_ysat(θsat,ϕ∞)
    za = PSP_za(ysat)
    ϕ0 = PSP_ϕ0(ϕ∞,za,ysat)
    asat = PSP_asat(xgs,G)
    ϕ2 = PSP_ϕ2(phit,G,Ds,asat)
    ϕsat = PSP_ϕsat(ϕ0,ϕ2)
    Vdsat = ϕsat-phit*log(1+(ϕsat*(ϕsat-2*asat*phit))/G^2/Ds/phit/phit)
    return Vdsat
end

function PSP_Vdseff(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox)
    Vdsat = PSP_Vdsat(Vg,Vs,Na,tox,Vfb,T,φb,Cox)
    AX = 3.0
    Vds = Vd-Vs
    return Vds/(1+(Vds/Vdsat)^AX)^(1/AX)
end

PSP_Gvsat(θsat,Δψ) = 0.5*(1+sqrt(1+2*(θsat*Δψ)^2))

PSP_T1(Vds,Vdse,Δψ,VP=0.05)=log(abs(1+(Vds-Δψ)/VP)/abs(1+(Vdse-Δψ)/VP))

PSP_ΔLL(T1,ALP=0.01) =ALP*T1

PSP_GΔL(ΔLL) = 1/(1+ΔLL+ΔLL^2)

function PSP_Gvsat(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox,CFD,CF)
    Vdseff = PSP_Vdseff(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox)
    Vds = Vd-Vs
    Vdsx = Vds^2/(sqrt(Vds^2+0.01)+0.1)
    ΔVgb = CF*2*(sqrt(1.0+CFD*Vdsx)-1.0)/CFD
    Vg = Vg+ΔVgb
    psiss,Es,Ds,Ps,xgs = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vs)
    psisd,Ed,Dd,Pd,xds = PSP_psiS(Vg,Na,tox,Vfb,T,φb,Vdseff)
    xds = (psisd-psiss)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    qis =  PSP_qis(G,Ds,xgs,Ps,phit)
    αs = PSP_αs(G,Es,Ps)
    qbs = PSP_qbs(phit,G,Ps)
    Eeff = PSP_Eeff(qbs,qis,Cox,ϵ₀*11.7)
    μe = PSP_μe(T)
    θμ = PSP_θμ(T)
    CS = 0
    Gmob = PSP_Gmob(μe,Eeff,θμ,CS,qbs,qis)
    wsat = PSP_wsat(qis)
    θsat = PSP_θsat(T,Gmob,wsat)
    ϕ∞ = PSP_ϕ∞(qis,PSP_αs(G,Es,Ps),phit)
    ysat = PSP_ysat(θsat,ϕ∞)
    za = PSP_za(ysat)
    ϕ0 = PSP_ϕ0(ϕ∞,za,ysat)
    asat = PSP_asat(xgs,G)
    ϕ2 = PSP_ϕ2(phit,G,Ds,asat)
    ϕsat = PSP_ϕsat(ϕ0,ϕ2)
    return Gmob*PSP_GΔL(PSP_ΔLL(PSP_T1(Vds,Vdseff,xds))),Vdseff,psiss,Es,Ds,Ps,xgs,psisd,Ed,Dd,Pd,xds,G,phit
end 

function PSP_ψms(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox,CFD,CF)
    Gmob,Vdseff,psiss,Es,Ds,Ps,xgs,psisd,Ed,Dd,Pd,xds,G,phit = PSP_Gvsat(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox,CFD,CF)
    φb *= 2.0
    x_s = psiss/phit
    x_d = psisd/phit
    xg = (Vg-Vfb)/phit
    xds = x_d-x_s
    if xg > 0.0
        xm = (x_s+x_d)*0.5
        Em = sqrt(Es*Ed)
        Dbar = (Ds+Dd)*0.5
        Dm = Dbar + xds^2/8.0*(Em-2/G^2)
        Pm = xm-1.0+Em
        xgm = G*sqrt(Dm+Pm)
        
        return xm*phit,Em,Dm,Pm,xds,xgm,Gmob
    else 
        return 0,0,0,0,0,0,0
    end
end

function PSP_Id(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox,CFD,CF,BETN=7e-2,U0=5e-2)
    β = BETN*U0*Cox*(T/300.0)^(-1.0)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    xm,Em,Dm,Pm,xds,xgm,Gmob = PSP_ψms(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox,CFD,CF)
    φb *= 2.0
    qim = PSP_qim(G,phit,Dm,Pm,xgm)
    α = PSP_α(G,Em,Pm)
    return β*(qim + phit*α)*phit*xds/Gmob
end

function PSP_Id(Vg,Vd,Vs,MOSf::MOSStructure,CFD,CF,BETN=7e-2,U0=5e-2)
   return PSP_Id(Vg,Vd,Vs,abs(C(MOSf)),MOSf.tox,VFB(MOSf),Temperature(MOSf),ϕb(MOSf),Cox(MOSf),CFD,CF,BETN,U0)
end

function PSP_ID(Vg,Vd,Vs,Na,tox,Vfb,T,φb,Cox)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    xm,Em,Dm,Pm,xds,xgm = PSP_ψm(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    φb *= 2.0
    qim = PSP_qim(G,phit,Dm,Pm,xgm)
    α = PSP_α(G,Em,Pm)
    return Cox*(qim + phit*α)*phit*xds
end

function PSP_ID(Vg,Vd,Vs,MOSf::MOSStructure)
   return PSP_ID(Vg,Vd,Vs,abs(C(MOSf)),MOSf.tox,VFB(MOSf),Temperature(MOSf),ϕb(MOSf),Cox(MOSf))
end

function qim_PSP(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    xm,Em,Dm,Pm,xds,xgm = PSP_ψm(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    return qim(G,phit,Dm,Pm,xgm)
end

function α_PSP(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    phit = kb*T
    G = γ(Na,tox)/sqrt(phit)
    xm,Em,Dm,Pm,xds,xgm = PSP_ψm(Vg,Vd,Vs,Na,tox,Vfb,T,φb)
    return α(G,Em,Pm)

end


