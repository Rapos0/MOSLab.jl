

struct MOSResult
    structure::MOSStructure
    result::Semiconductor1D
    VGB::Number
    ψs::Number
end

tox(MOS::MOSResult) = tox(MOS.structure)


C(MOS::MOSStructure) = C(MOS.semiconductor)
tox(MOS::MOSStructure) = MOS.tox
Temperature(MOS::MOSStructure) = Temperature(MOS.semiconductor)
VFB(MOS::MOSStructure) = φₘ(MOS.metal,Temperature(MOS)) - χ(MOS.semiconductor) - Ec(MOS.semiconductor)
ϕb(MOS::MOSStructure) = abs(Ef(MOS.semiconductor) - Ei(MOS.semiconductor))
Cox(MOS::MOSStructure) = MOS.oxide.ϵr*ϵ₀/tox(MOS)
n(MOS::MOSStructure,ψ=0.0) = n(MOS.semiconductor,ψ)
p(MOS::MOSStructure,ψ=0.0) = p(MOS.semiconductor,ψ)
ξ(MOS::MOSStructure,ψ=0.0) = ξ(MOS.semiconductor,ψ)
ξ(MOS::MOSStructure,ψ,ϕ) = ξ(MOS.semiconductor,ψ,ϕ)
Ei(MOS::MOSStructure) = Ei(MOS.semiconductor)
Ec(MOS::MOSStructure) = Ec(MOS.semiconductor)
Ev(MOS::MOSStructure) = Ev(MOS.semiconductor)
Nv(MOS::MOSStructure) = Nv(MOS.semiconductor)
Nc(MOS::MOSStructure) = Nc(MOS.semiconductor)
Ef(MOS::MOSStructure) = Ef(MOS.semiconductor)
ni(MOS::MOSStructure) = ni(MOS.semiconductor)
#Temperature(MOS::MOSStructure) = Temperature(MOS.semiconductor)
Rₗ(MOS::MOSStructure) = MOS.tsi
γ(MOS::MOSStructure) = γ(abs(C(MOS)),MOS.tox)
DistFun(MOS::MOSStructure,E) = DistFun(MOS.semiconductor,E)


function ψs₀(Vg,Na,tox,Vfb,T)
    Cox = ϵ₀*3.9/tox    
    ϕt = kb*T
    Gamma = Cox^-1*sqrt(2*q*11.7*ϵ₀*Na)
    α = -Gamma*sqrt(ϕt)
    ϕf = ϕt*log(Na/ni(T,BoltzmanDist(),Silicon()))
    ψs=0
    ug = Vg-Vfb
    if ug == 0.0
        return 0.0
    end
    uc = 2*ϕf
    ec = exp(-uc/ϕt)
    K = 1 - ec
    if uc > 0
        ϕsat = ug+Gamma^2/2.0*(1- sqrt(max(0,1+4*(ug-ϕt)/(Gamma^2))))
        if ϕsat >= uc && ug > 0
            ψs = uc + ϕt*log(max(1,(ug-uc)^2/Gamma^2/ϕt))
        else
            if ug < α
                ψs = -2*ϕt*log(ug/α)
            elseif ug > 0 && abs(ug) > ϕt 
                ψs = ϕsat
            end
        end
    end
    ψs = ug > 0 ? min(ug,ψs) : max(ug,ψs)
    return ψs
end


function Cap(Vg,MOS::MOSStructure)
    ψ = ψs(Vg,MOS)
    Xs = 1.0/Cs(ψ,MOS)
    return 1/(1/Cox(MOS)+Xs)
end

function Cs(ψ,MOS::MOSStructure)
    if ψ ∈ -0.03..0.03
        x0 = -0.03
        x1 = 0.03
        C0 = abs(Zygote.gradient(x->-ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,x),x0)[1])
        C1 = abs(Zygote.gradient(x->-ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,x),x1)[1])  
        α = (ψ-x0)/(x1-x0)
        return C0*(C1/C0)^α    #exp interp
    end  
    return abs(Zygote.gradient(x->-ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,x),ψ)[1])
end

function Qs(Vg,MOS::MOSStructure)
    return -ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,ψs(Vg,MOS))
end

ψs_PSP(Vg,MOS::MOSStructure) = PSP_VaricapPsi(Vg,abs(C(MOS)),MOS.tox,VFB(MOS),Temperature(MOS),ϕb(MOS))

ψs_PSPs(Vg,V,MOS::MOSStructure) = PSP_psiS(Vg,abs(C(MOS)),MOS.tox,VFB(MOS),Temperature(MOS),ϕb(MOS),V)[1]

function ψs(Vg,MOS::MOSStructure)
    ψ₀, = PSP_VaricapPsi(Vg,abs(C(MOS)),MOS.tox,VFB(MOS),Temperature(MOS),ϕb(MOS))
    QF(psi) = -ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,psi)
    fpsiB(psi) = psi + VFB(MOS) - QF(psi)/Cox(MOS)
    return nlsolve(x->fpsiB(x[1])-Vg,[ψ₀]).zero[1]
end
function ψs(Vg,ϕ,MOS::MOSStructure)
    ψ₀, = ψs_PSPs(Vg,ϕ,MOS)
    QF(psi) = -ϵ₀*ϵᵣ(MOS.semiconductor)*ξ(MOS.semiconductor,psi,ϕ)
    fpsiB(psi) = psi + VFB(MOS) - QF(psi)/Cox(MOS)
    return nlsolve(x->fpsiB(x[1])-Vg,[ψ₀]).zero[1]
end

function Poisson(Vg,MOS::MOSStructure,N=501)
    #println("New Version Build PN")
    Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(MOS.semiconductor,Rₗ(MOS))
    Ld = Ld0/lscal
    if Ld > 1.0
        println("Silicon Too small")
    end
    F(x) = DistFun(MOS,x)
    @variables x ψ(..)
    D = Differential(x)^2
    C0 = C(MOS)/Cscal
    ul = (Ef(MOS)-Ei(MOS)+Vg)/Vscal
    eqs = [λ*D(ψ(x)) ~ δ*F(ψ(x))-δ*F(-ψ(x))-C0]
    bcs = [ψ(0) ~ ψs(Vg,MOS)/Vscal+ul,
            ψ(Ld) ~ ul]
    domain = [x ∈ (0.0,Ld)]
    pde_system = PDESystem(eqs,bcs,domain,[x],[ψ(x)])
    discretization = MOLFiniteDifference([x=>Ld/(N-1)],nothing;centered_order=2)
    prob = discretize(pde_system,discretization)
    sol = solve(prob)*Vscal .+ Ei(MOS)
    x = collect(LinRange(0,Ld0,length(sol)))
    push!(x,Rₗ(MOS))
    push!(sol,0)  
    return MOSResult(MOS,Semiconductor1D(MOS.semiconductor,x,sol),Vg,sol[1])
end

function PoissonVoronoi(Vg,MOS::MOSStructure,N=20)
    Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(MOS.semiconductor,Rₗ(MOS))
    x,sol = PoissonVoronoi(Ec(MOS),Ev(MOS),Nc(MOS),Nv(MOS),0.0,ψs(Vg,MOS),0.0,lscal,Ld0,MOS.tsi,λ,N,N,x->C(MOS),Vscal,Cscal)
    return MOSResult(MOS,Semiconductor1D(MOS.semiconductor,x,sol),Vg,sol[1])
end

function BandDiagram(MOS::MOSResult)
    tsi = MOS.structure.tsi
    print(tsi)
    BDplot = BandDiagramJ(MOS.result,tsi)
    BDplot = plot!(BDplot,[-tox(MOS)*1e4,-tox(MOS)*1e4],[-10,10],label=nothing,c=:black)
    BDplot = plot!(BDplot,[0,0],[-10,10],label=nothing,c=:black,ylims=(-2.5,2.5))
    BDplot = plot!(BDplot,[tox(MOS)*1e4-0.1,-tox(MOS)*1e4],[-MOS.VGB+VFB(MOS.structure),-MOS.VGB+VFB(MOS.structure)],label=nothing,c=:black,xlims=(-2.0*tox(MOS)*1e4,tsi*1e4))
    BDplot = plot!(BDplot,[-tox(MOS)*1e4,0],[-MOS.VGB+VFB(MOS.structure),Ef(MOS.result[1])-MOS.ψs],label=L"\psi_{ox}",c=:red,xlabel=L"y \ [\mu m]",legend=:outertopright)
    return BDplot
end