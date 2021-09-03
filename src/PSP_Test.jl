

function sp_s(xg,Gf,phit,phib,Vsbstar)
    inv_phit1 = 1/phit
    oneSixth = 1/6.0
    invSqrt2 = 1/sqrt(2.0)
    xi           =  1.0 + Gf * invSqrt2
    inv_xi       =  1.0 /  xi
    Ux           =  Vsbstar * inv_phit1
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
end

P3(u) =  (1.0 + (u) * (1.0 + 0.5 * ((u) * (1.0 + (u)/3.0))))


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

function Vp_EKV(Vg,Vfb,Na,tox,T,phif,sigma =0.0 )
    GAMMA = γ(Na,tox)
    PHI = 2*phif
    VTO = Vfb+PHI+GAMMA*sqrt(PHI)
    VGprime = Vg - VTO + PHI + GAMMA * sqrt(PHI)
    VP = VGprime - PHI - GAMMA* (sqrt(VGprime+(GAMMA/2.0)*(GAMMA/2.0))-(GAMMA/2.0));
    return VP
end

function n_EKV(Vg,Vfb,Na,tox,T,phif,sigma =0.0 )
    GAMMA = γ(Na,tox)
    PHI = 2*phif
    VP = Vp_EKV(Vg,Vfb,Na,tox,T,phif)
    n = 1.0 + GAMMA / (2.0*sqrt(PHI + VP + 4.0*$vt));
    return n
end

function q_EKV(Vg,V,Vfb,Na,tox,T,phif,sigma =0.0 )
    VP = Vp_EKV(Vg,Vfb,Na,tox,T,phif)
    x=(VP-VS)/kb/T 
    iff = (ln(1.0+exp( x /2.0)))*(ln(1.0+exp( x /2.0)));
    return iff
end

function id_EKV(Vg,Vd,Vs,Vfb,Na,tox,T,phif,sigma =0.0 )
    pt = kb*T
    Ispec = 2*n_EKV(Vg,Vfb,Na,tox,T,phif)*pt*pt
    return Ispec*(q_EKV(Vg,Vs,Vfb,Na,tox,T,phif,sigma)-q_EKV(Vg,Vd,Vfb,Na,tox,T,phif,sigma))

end