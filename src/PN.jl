struct PNInput
    pInterval::Interval
    nInterval::Interval
    nMaterial::SemiconductorData
    pMaterial::SemiconductorData

end

PNInput(T,N,L) = PNInput(Domain(Interval{:close,:open}(0,L/2.0)),Domain(Interval{:close,:close}(L/2.0,L)),SemiconductorData(T,BoltzmanDist(),NSilicon(N,0)),SemiconductorData(T,BoltzmanDist(),PSilicon(N,0)))

function PoissonEq(PN::PNInput)
    F(x) = DistFun(PN.pMaterial,x)
    @variables x ψ(..)
    @parameters δ λ Cfun(x)
    Dxx = Differential(x)^2
    eqs = [
        λ*Dxx(ψ(x)) ~ δ*F(ψ(x))-p(ψ(x))-Cfun(x*lscal)    
    ]
    return eqs
end

function Poisson(Vg,PN::PNInput,Nres = 21.0)
    HeavisideC(x) = (1+sign(x+eps(Float64)))*0.5
    Heaviside(x) = (1+sign(x-eps(Float64)))*0.5
    Il(x,i::Interval) = HeavisideC(x-i.left)-Heaviside(x-i.right)
    Ir(x,i::Interval) = HeavisideC(x-i.left)-Heaviside(x-i.right)

    Vscal,lscalp,CN,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(PN.pMaterial,PN.pInterval.right)
    Vscal,lscaln,CP,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(PN.nMaterial,PN.nInterval.right)
    Cscal = max(abs.(CN),abs.(CP))
    lscal = max(lscaln,lscalp)
    f(x,a,b,Nn,Np) = C(Nn)*Il(x,a) + C(Np)*Ir(x,b) 
    g(x) = f(x,PN.pInterval,PN.nInterval,PN.pMaterial,PN.nMaterial)
    
    ul = -Ei(PN.pMaterial)/Vscal
    ur = -Ei(PN.nMaterial)/Vscal
    
    #C0 = Cdoping./Cscal
   
    dx = Ld0/Nres/lscal

    F(x) = DistFun(PN.pMaterial,x)
    n(ψ) = δ*F(ψ)
    p(ψ) = δ*F(-ψ)
    @variables x ψ(..)

    Dxx = Differential(x)^2
    Dx = Differential(x)


    Cfun(x) = g(x)/Cscal
    eqs = [
        λ*Dxx(ψ(x)) ~ n(ψ(x))-p(ψ(x))-Cfun(x*lscal)    
    ]
    bcs = [ψ(0) ~ ul+Vg/Vscal,
            ψ(1.0) ~ ur]
    domain = [x ∈ (0.0,1.0)]
    pde_system = PDESystem(eqs,bcs,domain,[x],[ψ(x)])
    discretization = MOLFiniteDifference([x=>dx],nothing;centered_order=2)
    prob = discretize(pde_system,discretization)
    L = length(prob.u0)
    L ÷ 2
    prob.u0[1:L÷2] .= ul
    prob.u0[1] = ul+ Vg/Vscal
    prob.u0[L÷2:L] .= ur
    prob.u0
    sol = solve(prob,jac=true)*Vscal
    if iseven(L)
        dxP = 0
    else
        dxP = dx*lscal
    end
    xP = PN.pInterval.left:dx*lscal:(PN.pInterval.right-dxP)
    LP = length(xP)
    SDP = Semiconductor1D(PN.pMaterial,xP,sol[1:LP].+Ei(PN.pMaterial))
    xN = PN.nInterval.left:dx*lscal:PN.nInterval.right
    LN = length(xN)
    SDN = Semiconductor1D(PN.nMaterial,xN,sol[LP+1:LP+LN].+Ei(PN.nMaterial))
    return Semiconductor1D(vcat(SDP.points,SDN.points))
    
end

function DD1D(Vg,PN::PNInput,Nres=21.0)
    HeavisideC(x) = (1+sign(x+eps(Float64)))*0.5
    Heaviside(x) = (1+sign(x-eps(Float64)))*0.5
    Il(x,i::Interval) = HeavisideC(x-i.left)-Heaviside(x-i.right)
    Ir(x,i::Interval) = HeavisideC(x-i.left)-Heaviside(x-i.right)

    Vscal,lscalp,CN,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(PN.pMaterial,PN.pInterval.right)
    Vscal,lscaln,CP,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(PN.nMaterial,PN.nInterval.right)
    Cscal = max(abs.(CN),abs.(CP))
    lscal = max(lscaln,lscalp)
    f(x,a,b,Nn,Np) = C(Nn)*Il(x,a) + C(Np)*Ir(x,b) 
    g(x) = f(x,PN.pInterval,PN.nInterval,PN.pMaterial,PN.nMaterial)
   
    τ = 1e-6
    τn = τ/tscal
    τp = τ/tscal
    
    un = 1400/μscal
    up = 450/μscal

    Rfunc(n,p) = δ/(τn*(p/δ+1)+τp*(n/δ+1))

    dx = Ld0/Nres/lscal

    F(x) = DistFun(PN.pMaterial,x)
    nf(ψ,ϕ) = δ*F(ψ-ϕ)
    pf(ψ,ϕ) = δ*F(-ψ+ϕ)
    
    @variables x u(..) ϕn(..) ϕp(..)

    Dxx = Differential(x)^2
    Dx = Differential(x)
    
    Cfun(x) = g(x)/Cscal
    ul = -Ei(PN.pMaterial)/Vscal
    ur = -Ei(PN.nMaterial)/Vscal
    
    #nD0 = (n(PN.pMaterial)+Nd(PN.pMaterial))/Cscal
    nD0 = (n(PN.pMaterial))/Cscal
    ψn0 = find_zero(x->nf(ul,x/Vscal)-nD0,(-100.0,10.0))*Vscal +  Vg/Vscal
    
    #pD0 = (p(PN.pMaterial)+Na(PN.pMaterial))/Cscal
    pD0 = (p(PN.pMaterial))/Cscal
    pD0 = ni(PN.pMaterial)
    ψp0 = find_zero(x->pf(ul,x/Vscal)-nD0,(-100.0,10.0))*Vscal + Vg/Vscal
    
    #nDL = (n(PN.nMaterial)+Nd(PN.nMaterial))/Cscal
    nDL = (n(PN.nMaterial))/Cscal
    ψnL = find_zero(x->nf(ul,x/Vscal)-nDL,(-100.0,100.0))*Vscal
    
    #pDL = (p(PN.nMaterial)+Na(PN.nMaterial))/Cscal
    pDL = (p(PN.nMaterial))/Cscal
    ψpL = find_zero(x->pf(ul,x/Vscal)-nDL,(-100.0,100.0))*Vscal
    
    
    eqs = [
        λ*Dxx(u(x)) ~ nf(u(x),0)-pf(u(x),0)-Cfun(x*lscal),
        Dxx(ϕn(x)) ~ -Rfunc(nf(u(x),ϕn(x)),pf(u(x),0.0))/nf(u(x),ϕn(x))/un+Dx(ϕn(x))*(Dx(ϕn(x)) - Dx(u(x))),
        Dxx(ϕp(x)) ~  Rfunc(nf(u(x),ϕn(x)),pf(u(x),ϕp(x)))/pf(u(x),ϕp(x))/up-Dx(ϕp(x))*(Dx(ϕp(x)) - Dx(u(x)))
        #Dx(up*p(ψ(x),ϕp(x)))*Dx(ϕp(x))+(up*p(ψ(x),ϕp(x)))*Dxx(ϕp(x)) ~ -R(n(ψ(x),ϕn(x)),p(ψ(x),ϕp(x)))
    ]

    bcs = [u(0.0) ~ ul,
           u(1.0) ~ ur,
           ϕn(0.0) ~ ψn0,
           ϕn(1.0) ~ ψnL,
           ϕp(0.0) ~ ψp0,
           ϕp(1.0) ~ ψpL,
            ]
        domain = [x ∈ (0.0,1.0)]
        pde_system = PDESystem(eqs,bcs,domain,[x],[u(x),ϕn(x),ϕp(x)])
        discretization = MOLFiniteDifference([x=>dx],nothing;centered_order=2)
        prob = discretize(pde_system,discretization)
        L = length(prob.u0)÷3
        prob.u0[1:L÷2] .= ul+ Vg/Vscal
        prob.u0[L÷2+1:L] .= ur
        sol = solve(prob,jac=true)*Vscal 
        if iseven(L)
            dxP = 0
        else
            dxP = dx*lscal
        end
        ψsol = sol[1:L]
        ϕnsol = sol[L+1:2L]
        ϕpsol = sol[2L+1:3L]

        xP = PN.pInterval.left:dx*lscal:(PN.pInterval.right-dxP)
        LP = length(xP)
        SDP = Semiconductor1D(PN.pMaterial,xP,ψsol[1:LP].+Ei(PN.pMaterial),ϕnsol[1:LP],ϕpsol[1:LP])
        xN = PN.nInterval.left:dx*lscal:PN.nInterval.right
        LN = length(xN)
        SDN = Semiconductor1D(PN.nMaterial,xN,ψsol[LP+1:LP+LN].+Ei(PN.nMaterial),ϕnsol[LP+1:LP+LN],ϕpsol[LP+1:LP+LN])
        
        return Semiconductor1D(vcat(SDP.points,SDN.points))
 
end