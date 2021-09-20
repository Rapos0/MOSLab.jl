struct PhysicalConfiguration
    Dom::Domain
    S::AbstractMaterial
end

struct PhysicalBoundaryConditions{T}
    position
    values::AbstractVector
    variable
end

struct SemiconductorProblem
    physicalConfiguration::Vector{PhysicalConfiguration}
    #bcs::AbstractVector{PhysicalBoundaryConditions{T}}
end

function DoppingFunction(dom::Domain,C::Number)
    return f(x) = (x ∈ dom)*C
end

# function RelaxStep(S::Semiconductor1D,x::AbstractVector)
#     applyPsiProfile(S,x)
#     return relaxationCostFunction(S)
# end

# function Relax(S::Semiconductor1D)
#     error(x) = RelaxStep(S,x)
#     return optimize(error,0.0*collect(1:length(S)),GradientDescent())
# end

function GummelPoisson(PP::Semiconductor1D,ψ₀,ψₗ,α=10.0;abserror=1e-12,itermax=100)
    n_i = ni(PP[1])
    ϕₜ = kb*Temperature(PP)
    Cdoping = C(PP)
    Cscal = maximum(vcat([n_i],abs.(Cdoping)))
    lscal = Rₗ(PP)
    Nnodes = length(PP)
    ui = -Ei(PP)/ϕₜ
    return GummelPoisson(ψ₀/ϕₜ + ui[1],ψₗ/ϕₜ + ui[end],ui,Cdoping,Nnodes,n_i,lscal,ϕₜ,r(PP),α;abserror=abserror,itermax=itermax)*ϕₜ
    return r0*lscal,u*ϕₜ,e
end

function GummelPoisson(u₀,uₗ,ui,C,Nnodes,n_i,lscal,Vscal,r0,α=1.0;abserror=1e-12,itermax=100,relerror=1e-12)    
    Cscal = maximum(vcat([n_i],abs.(C)))
    λ = ϵ₀*11.7*Vscal/(lscal^2*q*Cscal)
    r0 ./=lscal
    C = C./Cscal
    δ = n_i/Cscal
    dr = diff(r0)
    #Operators 
    ∇² = CenteredDifference(2,5,dr,Nnodes-2)
    Bc =  DirichletBC(u₀,uₗ)
    φBC = DirichletBC(0.0,0.0)
    ## Initialization
    u0 = ui
    u0[1] = u₀
    u0[end] = uₗ
    ## Step
    u1 = copy(u0)
    limexp(x) = x > 1e2 ? exp(1e2) : exp(x)
    n(u) = δ*limexp(u)
    p(u) = δ*limexp(-u)
    R(u) = -λ*(∇²*Bc*u)+δ*limexp.(u)-δ*limexp.(-u) - C[2:end-1]
    e = Inf64
    iter = 1
    Rmin = Inf64
    uf = similar(u0)
    nR = Inf64
    while e > abserror && iter < itermax && nR > relerror
        res = R(u0[2:end-1])
        SGL(u) = -λ*(∇²*φBC*u) + (n.(u)+p.(u)).*u + res
        ϕ = zeros(length(u0))
        ϕ[2:end-1] .= nlsolve(SGL,ϕ[2:end-1],autodiff = :forward).zero
        b = abs.(ϕ) .> 1.0
        # This damping procedure is inspired from Solid-State Electronics, vol. 19,
        # pp. 991-992 (1976).
        ϕ[b] = log.(1.0.+abs.(ϕ[b]).*1.72).*sign.(ϕ[b])
        u1[2:end-1] .+= ϕ[2:end-1] 
        nR = sum(abs.(ϕ[2:end-1]))/length(res)
        e = norm(res)
        #println("i:: $(iter) norm:$(nR) res:$(e)")
        if(nR < Rmin)
            uf = u0
            Rmin = nR
        end
        u0 = u1
        iter = iter+1
        
    end
    return uf
end

function scallingNv(S::SemiconductorData,L)
    n_i = ni(S)
    ϕₜ = kb*Temperature(S)
    Vscal = ϕₜ 
    Cscal = max(n_i,abs.(C(S)),Nv(S))
    lscal = L
    λ = ϵ₀*11.7*Vscal/(lscal^2*q*Cscal)
    δ = n_i/Cscal
    μscal = 1000
    Ld0 = sqrt(ϵ₀*11.7* ϕₜ/q/abs(C(S)))
    Rscal = ϕₜ*μscal*Cscal/lscal/lscal
    tscal = lscal^2/ϕₜ/μscal
    Jscal = q*Vscal*μscal*Cscal/lscal
    Dscal = μscal*Vscal
    return Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i

end

function scallingLD(S::SemiconductorData,L)
    n_i = ni(S)
    ϕₜ = kb*Temperature(S)
    Vscal = ϕₜ 
    Cscal = max(n_i,abs.(C(S)))
    lscal = sqrt(ϵ₀*11.7* ϕₜ/q/abs(C(S)))
    λ = ϵ₀*11.7*Vscal/(lscal^2*q*Cscal)
    δ = n_i/Cscal
    μscal = 1000
    Ld0 = sqrt(ϵ₀*11.7* ϕₜ/q/abs(C(S)))
    Rscal = ϕₜ*μscal*Cscal/lscal/lscal
    tscal = lscal^2/ϕₜ/μscal
    Jscal = q*Vscal*μscal*Cscal/lscal
    Dscal = μscal*Vscal
    return Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i

end

function scallingNvLD(S::SemiconductorData,L)
    n_i = ni(S)
    ϕₜ = kb*Temperature(S)
    Vscal = ϕₜ 
    Cscal = max(n_i,abs.(C(S)),Nv(S))
    lscal = sqrt(ϵ₀*11.7* ϕₜ/q/abs(C(S)))
    λ = ϵ₀*11.7*Vscal/(lscal^2*q*Cscal)
    δ = n_i/Cscal
    μscal = 1000
    Ld0 = sqrt(ϵ₀*11.7* ϕₜ/q/abs(C(S)))
    Rscal = ϕₜ*μscal*Cscal/lscal/lscal
    tscal = lscal^2/ϕₜ/μscal
    Jscal = q*Vscal*μscal*Cscal/lscal
    Dscal = μscal*Vscal
    return Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i

end

function limexp(x,lm=-50,lM=50)
    if x < lm
        return exp(lm)*((x-lm)+1.0)
    elseif x > lM
        return exp(lM)*((x-lM)+1.0)
    else
        return exp(x)
    end

end
function DDVoronoiMOSFET(X,Ec,Ev,Nc,Nv,Vbi,Vg,Vd,lscal,Rscal,Cscal,Vscal,λ,Jn,Jp,R,C,fetIn::MOSFETInputDeck;u0=nothing,n0=nothing,p0=nothing,ξ₀=1e-5,verbose=false)
    ## Selective Space Discretization
    
    #Lint = maximum(X)
    ## Boltzman Dist
    F(x) = limexp(x)
    ηₙ(u,un) = u-Ec/Vscal
    ncf(u,un) = Nc*F(ηₙ(u,un))/Cscal
    ηₚ(u,up) = -u+Ev/Vscal
    pcf(u,up) = Nv*F(ηₚ(u,up))/Cscal
    nc(u,un) = un
    pc(u,up) = up

    #Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])/meas(edge)
                
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))/meas(edge)
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))/meas(edge)
                # @info "u $(map(x->x.value,u.val)) - m: $(meas(edge)) - f:$(map(x->x.value,f.val))"
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal,node[2]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(X,physics)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])
    
    Vs = 0.0
    us = 0.0
    ug = (Vbi+Vg)/Vscal
    ud = (Vd)/Vscal
    uns = ncf(us,0.0)*exp(-Vs/Vscal)  
    ung = ncf(ug,0.0)*exp(-Vg/Vscal)  
    und = ncf(ud,0.0)exp(-Vd/Vscal)
    ups = pcf(us,0.0)*exp(Vs/Vscal)
    upg = pcf(ug,0.0)exp(Vg/Vscal)
    upd = pcf(ud,0.0)exp(Vd/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,1,us)
    boundary_dirichlet!(sys,1,6,us)
    boundary_dirichlet!(sys,1,7,ug)
    boundary_dirichlet!(sys,1,8,ud)
  
    boundary_dirichlet!(sys,2,1,uns)
    boundary_dirichlet!(sys,2,6,uns)
    boundary_dirichlet!(sys,2,7,ung)
    boundary_dirichlet!(sys,2,8,und)

    boundary_dirichlet!(sys,3,1,ups)
    boundary_dirichlet!(sys,3,6,ups)
    boundary_dirichlet!(sys,3,7,upg)
    boundary_dirichlet!(sys,3,8,upd)

    uv = [us,ug,ud]
    nv = [uns,ung,und]
    pv = [ups,upg,upd]
    Cords = X[Coordinates]
    # Create a solution array
    inival=unknowns(sys,inival=0.1)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        for i in 1:size(inival)[2]
            inival[1,i] = dot(belongVector(Cords[1,i]*lscal,Cords[2,i]*lscal,fetIn),uv)
        end
    else
        inival[1,:] .= u0
    end
    if isnothing(n0)
        for i in 1:size(inival)[2]
            inival[2,i] = dot(belongVector(Cords[1,i]*lscal,Cords[2,i]*lscal,fetIn),nv)
        end
    else
        inival[2,:] = n0
    end
    if isnothing(p0)
        for i in 1:size(inival)[2]
            inival[3,i] = dot(belongVector(Cords[1,i]*lscal,Cords[2,i]*lscal,fetIn),pv)
        end
    else
        inival[3,:] .=p0 
    end
    
    solution=unknowns(sys)

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    #control.max_lureuse = 0
    control.verbose=verbose
    control.max_iterations = 500
    #control.Δu_opt=100
    control.damp_initial= ξ₀
    # Stationary solution of the problem
    VoronoiFVM.solve!(solution,inival,sys, control=control)
    
    return solution[1,:]*Vscal,solution[2,:]*Cscal,solution[3,:]*Cscal
end


function DDVoronoiMOSFET(X,Ec,Ev,Nc,Nv,Vbi,Vg,Vd,lscal,Rscal,Cscal,Vscal,Ld0,L,λ,nld,Jn,Jp,nn,R,C;u0=nothing,n0=nothing,p0=nothing,ξ₀=1e-5,verbose=false)
    ## Selective Space Discretization
    
    #Lint = maximum(X)
    ## Boltzman Dist
    F(x) = abs(x) > 100 ? exp(sign(x)*100) : exp(x)
    ηₙ(u,un) = u-Ec/Vscal
    ncf(u,un) = Nc*F(ηₙ(u,un))/Cscal
    ηₚ(u,up) = -u+Ev/Vscal
    pcf(u,up) = Nv*F(ηₚ(u,up))/Cscal
    nc(u,un) = un
    pc(u,up) = up

    #Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal,node[2]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(X,physics,unknown_storage=:sparse)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])
    
    Vs = 0.0
    us = 0.0
    ug = (Vbi+Vg)/Vscal
    ud = (Vd)/Vscal
    uns = ncf(us,0.0)*exp(-Vs/Vscal)  
    ung = ncf(ug,0.0)*exp(-Vg/Vscal)  
    und = ncf(ud,0.0)exp(-Vd/Vscal)
    ups = pcf(us,0.0)*exp(Vs/Vscal)
    upg = pcf(ug,0.0)exp(Vg/Vscal)
    upd = pcf(ud,0.0)exp(Vd/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,1,us)
    boundary_dirichlet!(sys,1,6,us)
    boundary_dirichlet!(sys,1,7,ug)
    boundary_dirichlet!(sys,1,8,ud)
  
    boundary_dirichlet!(sys,2,1,uns)
    boundary_dirichlet!(sys,2,6,uns)
    boundary_dirichlet!(sys,2,7,ung)
    boundary_dirichlet!(sys,2,8,und)

    boundary_dirichlet!(sys,3,1,ups)
    boundary_dirichlet!(sys,3,6,ups)
    boundary_dirichlet!(sys,3,7,upg)
    boundary_dirichlet!(sys,3,8,upd)


    # Create a solution array
    inival=unknowns(sys,inival=0.1)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        inival[1,1:L] .= us
        inival[1,L+1:2L] .= ud
    else
        inival[1,:] .= u0
    end
    if isnothing(n0)
        inival[2,1:L] .= uns
        inival[2,L+1:2L] .= und
        
    else
        inival[2,:] .= n0
    end
    if isnothing(p0)
        inival[3,1:L] .= ups
        inival[3,L+1:2L] .= upd
    else
        inival[3,:] .= p0
    end
    
    solution=unknowns(sys)

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    #control.max_lureuse = 0
    control.verbose=verbose
    control.max_iterations = 500
    #control.Δu_opt=100
    control.damp_initial= ξ₀
    # Stationary solution of the problem
    VoronoiFVM.solve!(solution,inival,sys, control=control)
    
    return solution[1,:]*Vscal,solution[2,:]*Cscal,solution[3,:]*Cscal
end

function DDVoronoiHC2D(X,Ec,Ev,Nc,Nv,Vbi,Va,Vb,lscal,Rscal,Cscal,Vscal,Ld0,L,λ,nld,Jn,Jp,nn,R,C,u0=nothing,n0=nothing,p0=nothing)
    ## Selective Space Discretization
    
    #Lint = maximum(X)
    ## Boltzman Dist
    F(x) = abs(x) > 100 ? exp(sign(x)*100) : exp(x)
    ηₙ(u,un) = u-Ec/Vscal
    ncf(u,un) = Nc*F(ηₙ(u,un))/Cscal
    ηₚ(u,up) = -u+Ev/Vscal
    pcf(u,up) = Nv*F(ηₚ(u,up))/Cscal
    nc(u,un) = un
    pc(u,up) = up

    #Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])/meas(edge)
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))/meas(edge)
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))/meas(edge)
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(X,physics,unknown_storage=:sparse)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])

    ul = 0.0+Va/Vscal
    ur = (Vbi+Vb)/Vscal
    unr = ncf(ul,0.0)*exp(-Va/Vscal)  
    unl = ncf(ur,0.0)exp(-Vb/Vscal)
    upr = pcf(ul,0.0)*exp(Va/Vscal)
    upl = pcf(ur,0.0)exp(Vb/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,4,ul)
    boundary_dirichlet!(sys,1,2,ur)
    boundary_dirichlet!(sys,2,4,unr)
    boundary_dirichlet!(sys,2,2,unl)
    boundary_dirichlet!(sys,3,4,upr)
    boundary_dirichlet!(sys,3,2,upl)

    # Create a solution array
    inival=unknowns(sys,inival=0.1)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        inival[1,1:L] .= ul
        inival[1,L+1:2L] .= ur
    else
        inival[1,:] .= u0
    end
    if isnothing(n0)
        inival[2,:] .= n0
    else
        inival[2,1:L] .= unr
        inival[2,L+1:2L] .= unl
    end
    if isnothing(p0)
        inival[3,1:L] .= upr
        inival[3,L+1:2L] .= upl
    else
        
    end
    
    solution=unknowns(sys)

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    #control.max_lureuse = 0
    control.verbose=true
    control.max_iterations = 500
    #control.Δu_opt=100
    control.damp_initial=1e-10
    # Stationary solution of the problem
    VoronoiFVM.solve!(solution,inival,sys, control=control)
    
    return solution[1,:]*Vscal,solution[2,:],solution[3,:]
end


function DDVoronoiHC(Ec,Ev,Nc,Nv,Vbi,Va,Vb,lscal,Rscal,Cscal,Vscal,Ld0,L,λ,nld,Jn,Jp,nn,R,C,u0=nothing)
    ## Selective Space Discretization
    dx0 = 3*Ld0/L
    dxld = dx0/nld
    dxS0 = dx0/nn
    p0 = 0
    p1 = dx0
    p2 = 0.5-dx0
    p3 = 0.5+dx0
    p4 = 1.0-dx0
    p5 = 1.0
    X1 = p0:dxld:p1
    X2 = p1:dxS0:p2
    X3 = p2:dxld:p3
    X4 = p3:dxS0:p4
    X5 = p4:dxld:p5
    X = vcat(X1,X2[2:end],X3[2:end],X4[2:end],X5[2:end])*L/lscal
    Lint = maximum(X)
    ## Boltzman Dist
    F(x) = abs(x) > 100 ? exp(sign(x)*100) : exp(x)
    ηₙ(u,un) = u-Ec/Vscal
    ncf(u,un) = Nc*F(ηₙ(u,un))/Cscal
    ηₚ(u,up) = -u+Ev/Vscal
    pcf(u,up) = Nv*F(ηₚ(u,up))/Cscal
    nc(u,un) = un
    pc(u,up) = up

    Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])/meas(edge)
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))/meas(edge)
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))/meas(edge)
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(Gx,physics,unknown_storage=:sparse)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])

    ul = 0.0+Va/Vscal
    ur = (Vbi+Vb)/Vscal
    unr = ncf(ul,0.0)*exp(-Va/Vscal)  
    unl = ncf(ur,0.0)exp(-Vb/Vscal)
    upr = pcf(ul,0.0)*exp(Va/Vscal)
    upl = pcf(ur,0.0)exp(Vb/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,1,ul)
    boundary_dirichlet!(sys,1,2,ur)
    boundary_dirichlet!(sys,2,1,unr)
    boundary_dirichlet!(sys,2,2,unl)
    boundary_dirichlet!(sys,3,1,upr)
    boundary_dirichlet!(sys,3,2,upl)

    # Create a solution array
    inival=unknowns(sys,inival=0.1)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        inival[1,1:L] .= ul
        inival[1,L+1:2L] .= ur
    else
        inival[1,:] .= u0
    end
    inival[2,1:L] .= unr
    inival[2,L+1:2L] .= unl
    inival[3,1:L] .= upr
    inival[3,L+1:2L] .= upl
    
    solution=unknowns(sys)

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    #control.max_lureuse = 0
    control.verbose=true
    control.max_iterations = 500
    #control.Δu_opt=100
    control.damp_initial=1e-10
    # Stationary solution of the problem
    VoronoiFVM.solve!(solution,inival,sys, control=control)
    
    return X*lscal,solution[1,:]*Vscal,solution[2,:],solution[3,:]
end


function DDVoronoiH(Ec,Ev,Nc,Nv,Vbi,Va,Vb,lscal,Ld0,L,λ,nld,Jn,Jp,nn,R,C,u0=nothing)
    ## Selective Space Discretization
    dx0 = 3*Ld0/L
    dxld = dx0/nld
    dxS0 = dx0/nn
    p0 = 0
    p1 = dx0
    p2 = 0.5-dx0
    p3 = 0.5+dx0
    p4 = 1.0-dx0
    p5 = 1.0
    X1 = p0:dxld:p1
    X2 = p1:dxS0:p2
    X3 = p2:dxld:p3
    X4 = p3:dxS0:p4
    X5 = p4:dxld:p5
    X = vcat(X1,X2[2:end],X3[2:end],X4[2:end],X5[2:end])*L/lscal
    Lint = maximum(X)
    ## Boltzman Dist
    F(x) = abs(x) > 100 ? exp(sign(x)*100) : exp(x)
    ηₙ(u,un) = u-Ec/Vscal
    nc(u,un) = Nc*F(ηₙ(u,un))*un/Cscal
    ηₚ(u,up) = -u+Ev/Vscal
    pc(u,up) = Nv*F(ηₚ(u,up))*up/Cscal

    Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])/meas(edge)
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))/meas(edge)
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))/meas(edge)
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(Gx,physics,unknown_storage=:sparse)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])

    ul = 0.0+Va/Vscal
    ur = (Vbi+Vb)/Vscal
    unr = exp(-Va/Vscal)  
    unl = exp(-Vb/Vscal)
    upr = exp(Va/Vscal)
    upl = exp(Vb/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,1,ul)
    boundary_dirichlet!(sys,1,2,ur)
    boundary_dirichlet!(sys,2,1,unr)
    boundary_dirichlet!(sys,2,2,unl)
    boundary_dirichlet!(sys,3,1,upr)
    boundary_dirichlet!(sys,3,2,upl)

    # Create a solution array
    inival=unknowns(sys,inival=0.1)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        inival[1,1:L] .= ul
        inival[1,L+1:2L] .= ur
    else
        inival[1,:] .= u0
    end
    inival[2,1:L] .= unr
    inival[2,L+1:2L] .= unl
    inival[3,1:L] .= upr
    inival[3,L+1:2L] .= upl
    
    solution=unknowns(sys)

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    #control.max_lureuse = 0
    control.verbose=true
    control.max_iterations = 500
    #control.Δu_opt=100
    control.damp_initial=1e-10
    # Stationary solution of the problem
    solve!(solution,inival,sys, control=control)
    
    return X*lscal,solution[1,:]*Vscal,-log.(solution[2,:])*Vscal,log.(solution[3,:])*Vscal
end

function DDVoronoi(Ec,Ev,Nc,Nv,Vbi,Va,Vb,lscal,Ld0,L,λ,nld,Jn,Jp,nn,R,C,control,u0=nothing)
    ## Selective Space Discretization
    dx0 = 3*Ld0/L
    dxld = dx0/nld
    dxS0 = dx0/nn
    p0 = 0
    p1 = dx0
    p2 = 0.5-dx0
    p3 = 0.5+dx0
    p4 = 1.0-dx0
    p5 = 1.0
    X1 = p0:dxld:p1
    X2 = p1:dxS0:p2
    X3 = p2:dxld:p3
    X4 = p3:dxS0:p4
    X5 = p4:dxld:p5
    X = vcat(X1,X2[2:end],X3[2:end],X4[2:end],X5[2:end])*L/lscal
    Lint = maximum(X)
    ## Boltzman Dist
    F(x) = abs(x) > 100 ? exp(sign(x)*100) : exp(x)
    ηₙ(u,un) = u-un-Ec/Vscal
    nc(u,un) = Nc*F(ηₙ(u,un))/Cscal
    ηₚ(u,up) = up-u+Ev/Vscal
    pc(u,up) = Nv*F(ηₚ(u,up))/Cscal

    Gx = VoronoiFVM.Grid(X)

    #Create PhysicalModel
   
    
 
    physics=VoronoiFVM.Physics(
        flux= function flux!(f,u,edge)
                f[1]=λ*(u[1,1]-u[1,2])/meas(edge)
                l = 2
                k = 1
                bm = fbernoulli(-(u[1,l]-u[1,k]))
                bp = fbernoulli((u[1,l]-u[1,k]))
                f[2]= -Jn*(bm*nc(u[1,k],u[2,k])-bp*nc(u[1,l],u[2,l]))/meas(edge)
                f[3]= Jp*(bp*pc(u[1,k],u[3,k])-bm*pc(u[1,l],u[3,l]))/meas(edge)
                #bpp,bmp=fbernoulli_pm(u[1,1]-u[1,-1])
                #f[3]= 0.0#Jp*(bp*pc(u[1,1],u[3,1])-bm*pc(u[1,2],u[3,2]))
        end,
        source=function (f,node)
            f[1] = C(node[1]*lscal)/Cscal
            f[2] = 0.0
            f[3] = 0.0
        end,
        reaction=function reaction!(f,u,node)
            f[1] = nc(u[1],u[2])-pc(u[1],u[3])
            f[2] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
            f[3] = R(nc(u[1],u[2])*Cscal,pc(u[1],u[3])*Cscal)/Rscal
        end)

    sys=VoronoiFVM.System(Gx,physics,unknown_storage=:sparse)

    # Add species 1 to region 1
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_species!(sys,3,[1])

    ul = 0.0+Va/Vscal
    ur = (Vbi+Vb)/Vscal
    unr = exp(-Va/Vscal)  
    unl = exp(-Vb/Vscal)
    upr = exp(Va/Vscal)
    upl = exp(Vb/Vscal)
    # Set boundary conditions
    boundary_dirichlet!(sys,1,1,ul)
    boundary_dirichlet!(sys,1,2,ur)
    boundary_dirichlet!(sys,2,1,Va/Vscal)
    boundary_dirichlet!(sys,2,2,Vb/Vscal)
    boundary_dirichlet!(sys,3,1,Va/Vscal)
    boundary_dirichlet!(sys,3,2,Vb/Vscal)


    # Create a solution array
    inival=unknowns(sys,inival=0.001)
    print(size(inival))
    L = length(inival[1,:])÷2
    if isnothing(u0)
        inival[1,1:L] .= ul
        inival[1,L+1:2L] .= ur
    else
        inival[1,:] .= u0
    end

    
    solution=unknowns(sys)

    # Create solver control info
   
    # Stationary solution of the problem
    VoronoiFVM.solve!(solution,inival,sys, control=control)
    
    return X*lscal,solution[1,:]*Vscal,-log.(solution[2,:])*Vscal,log.(solution[3,:])*Vscal
end




function PoissonVoronoi(Ec,Ev,Nc,Nv,Vbi,Va,Vb,lscal,Ld0,L,λ,nld,nn,C,Vscal,Cscal)
    ## Selective Space Discretization
   ## Selective Space Discretization
   dx0 = 3*Ld0/L
   dxld = dx0/nld
   dxS0 = dx0/nn
   p0 = 0
   p1 = dx0
   p2 = 0.5-dx0
   p3 = 0.5+dx0
   p4 = 1.0-dx0
   p5 = 1.0
   X1 = p0:dxld:p1
   X2 = p1:dxS0:p2
   X3 = p2:dxld:p3
   X4 = p3:dxS0:p4
   X5 = p4:dxld:p5
   X = vcat(X1,X2[2:end],X3[2:end],X4[2:end],X5[2:end])*L/lscal
   Lint = maximum(X)
   ## Boltzman Dist
   F(x) = abs(x) > 500 ? exp(sign(x)*500) : exp(x)
   ηₙ(u) = u-Ec/Vscal
   nc(u) = Nc*F(ηₙ(u))/Cscal
   ηₚ(u) = -u+Ev/Vscal
   pc(u) = Nv*F(ηₚ(u))/Cscal

   Gx = VoronoiFVM.Grid(X)

   #Create PhysicalModel
  
   

   physics=VoronoiFVM.Physics(
       flux= function flux!(f,u,edge)
               f[1]=λ*(u[1,1]-u[1,2])
       end,
       source=function (f,node)
           f[1] = C(node[1]*lscal)/Cscal
       end,
       reaction=function reaction!(f,u,node)
           f[1] = nc(u[1])-pc(u[1])
       end)

   sys=VoronoiFVM.System(Gx,physics,unknown_storage=:sparse)

   # Add species 1 to region 1
   enable_species!(sys,1,[1])


   ul = 0.0+Va/Vscal
   ur = (Vbi+Vb)/Vscal

   # Set boundary conditions
   boundary_dirichlet!(sys,1,1,ul)
   boundary_dirichlet!(sys,1,2,ur)
   
   # Create a solution array
   inival=unknowns(sys,inival=1.0)
   print(size(inival))
   L = length(inival[1,:])÷2
   inival[1,1:L] .= ul
   inival[1,L+1:2L] .= ur
   
   solution=unknowns(sys)

   # Create solver control info
   control=VoronoiFVM.NewtonControl()
   #control.max_lureuse = 0
   control.verbose=true
   control.max_iterations = 500
   #control.Δu_opt=100
   control.damp_initial=1e-2
   control.tol_absolute = 1e-4
   control.tol_relative = 1e-4
   control.damp_growth = 2.0
   # Stationary solution of the problem
   VoronoiFVM.solve!(solution,inival,sys, control=control)
   
    return X*lscal,solution[1,:]*Vscal
end
