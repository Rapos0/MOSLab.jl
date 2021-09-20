"""
    Semiconductor1D

    Represents a series on SemiconductorData objects with positions 
"""
mutable struct Semiconductor1D{K} <: AbstractArray{SemiconductorPoint{K},1} 
    points::AbstractArray{SemiconductorPoint{K},1}
end

Base.size(S::Semiconductor1D) = size(S.points)
Base.push!(S::Semiconductor1D,e::SemiconductorPoint{Float64}) = (S.points = vcat(S.points,[e]))
Base.length(S::Semiconductor1D) = length(S.points)
Base.getindex(S::Semiconductor1D,i) = S.points[i]
Base.lastindex(S::Semiconductor1D) = length(S)
Base.setindex!(S::Semiconductor1D,v,i) = (S.points[i] = v)

#Semiconductor1D(V::AbstractArray )= Semiconductor1D(V)
"""
    Semiconductor1D(S::SemiconductorData,r,psi) 

    Create a Semiconductor1D with S as the base semiconductor at positions r and electrostatic potential psi (referenced at the Fermi Level) at positions r (in cm)

    ```jldoctest
    julia> 
        T = 300 # K
        N_a = 5e17 # cm^-3
        r = LinRange(0,1e-4,101) # in cm
        psi = LinRange(0.0,1.0,101) # in V
        semi1D = Semiconductor1D(SemicondorData(T,BoltzmanDist(),PSilicon(N_a)),r,psi)
    ```
"""
Semiconductor1D(S::SemiconductorData,r,psi) = Semiconductor1D(SemiconductorPoint1D.(Ref(S),r,psi))

"""
    Semiconductor1D(S::SemiconductorData,r,psi,ϕn,ϕp) 

    Create a Semiconductor1D with S as the base semiconductor at positions r and potentials psi, electron quasi-Fermi level ϕn, and holes quasi-Fermi level ϕn (referenced at the Fermi Level) at positions r
"""
Semiconductor1D(S::SemiconductorData,r,psi,ϕn,ϕp) = Semiconductor1D(SemiconductorPoint1D.(Ref(S),r,psi,ϕn,ϕp))


Semiconductor1D(S::SemiconductorData,r::AbstractVector{T}) where T <: Number = Semiconductor1D(SemiconductorPoint1D.(Ref(S),r,zeros(length(r))))

#Semiconductor1D() = Semiconductor1D(Vector{SemiconductorPoint{Float64}})
"""
    Eg(S::Semiconductor1D) Bandgap at different positions
"""
Eg(S::Semiconductor1D) = Eg.(S.points)
"""
    Ec(S::Semiconductor1D) Conduction Band at different positions referenced to the Fermi Level
"""
Ec(S::Semiconductor1D) = Ec.(S.points)

"""
    EV(S::Semiconductor1D) Valence Band at different positions referenced to the Fermi Level
"""
Ev(S::Semiconductor1D) = Ev.(S.points)
Ef(S::Semiconductor1D) = Ef.(S.points)
"""
    Ei(S::Semiconductor1D) Intrinsic Fermi Level at different positions referenced to the Fermi Level
"""
Ei(S::Semiconductor1D) = Ei.(S.points)
"""
    n(S::Semiconductor1D) Electrons concentration at different positions referenced to the Fermi Level
"""
n(S::Semiconductor1D) = n.(S.points)
"""
    p(S::Semiconductor1D) Holes concentration at different positions in cm^-3
"""
p(S::Semiconductor1D) = p.(S.points)
n(S::Semiconductor1D,ψ) = n.(S.points,ψ)
p(S::Semiconductor1D,ψ) = p.(S.points,ψ)
"""
    Na(S::Semiconductor1D) Acceptors concentration at different positions in cm^-3
"""
Na(S::Semiconductor1D) = Na.(S.points)
"""
    Nd(S::Semiconductor1D) Donors concentration at different positions in cm^-3
"""
Nd(S::Semiconductor1D) = Nd.(S.points)
N(S::Semiconductor1D) = N.(S.points)
N(S::Semiconductor1D,ψ) = N.(S.points,ψ)
P(S::Semiconductor1D,ψ) = P.(S.points,ψ)
P(S::Semiconductor1D) = P.(S.points)
"""
    ni(S::Semiconductor1D) Intrinsic carriers concentration at different positions in cm^-3
"""
ni(S::Semiconductor1D) = ni.(S.points)
"""
    φf(S::Semiconductor1D) Ef-Ei at different positions in V
"""
φf(S::Semiconductor1D) = φf.(S.points)
"""
    r(S::Semiconductor1D) Positions in which the 1D Semiconductor is defined in cm
"""
r(S::Semiconductor1D) = r.(S.points)
"""
    Rₗ(S::Semiconductor1D) End positions in which the 1D Semiconductor is defined in cm
"""
Rₗ(S::Semiconductor1D) = maximum(r(S))
"""
    R₀(S::Semiconductor1D) Begin positions in which the 1D Semiconductor is defined in cm
"""
R₀(S::Semiconductor1D) = minimum(r(S))
"""
    L(S::Semiconductor1D) Semiconductor1D length in cm 
"""
L(S::Semiconductor1D) = Rₗ(S)-R₀(S)
"""
    ρ(S::Semiconductor1D) Semiconductor1D electric charge density
"""
ρ(S::Semiconductor1D) = ρ.(S.points)
#ρ(S::Semiconductor1D,ψ) = -q/ϵ₀/ϵᵣ(S[1]) *  (P(S,ψ).-N(S,ψ)) *S.W*S.W
"""
    ξ(S::Semiconductor1D) Semiconductor1D electric field calculated analytically
"""
ξ(S::Semiconductor1D) = ξ.(S.points)
ξ(S::Semiconductor1D,ψ) = ξ.(S.points,ψ)
ξ(S::Semiconductor1D,ψ,ϕ) = ξ.(S.points,ψ,ϕ)
"""
    ψ(S::Semiconductor1D) Electrostatic potential at the different positions
"""
ψ(S::Semiconductor1D) = ψ.(S.points)
Δχ(S::Semiconductor1D,δχ) = Δχ.(S.points,δχ)
"""
    Temperature(S::Semiconductor1D) Latice Temperature
"""
Temperature(S::Semiconductor1D) = Temperature(S.points[1])
"""
    ϵᵣ(S::Semiconductor1D) Relative permitivity
"""
ϵᵣ(S::Semiconductor1D)=ϵᵣ(S.points[1])
"""
    χ(S::Semiconductor1D) Work function
"""
χ(S::Semiconductor1D) = χ.(S.points)
"""
    Nv(S::Semiconductor1D) Effective density of states on the Valence Band
"""
Nv(S::Semiconductor1D) = Nv(S.points[1])
"""
    Nc(S::Semiconductor1D) Effective density of states on the Conduction Band
"""
Nc(S::Semiconductor1D) = Nc(S.points[1])
Pt(S::Semiconductor1D) = Pt.(S.points)
Nt(S::Semiconductor1D) = Nt.(S.points)
"""
    DistFun(S::Semiconductor1D) Distribution Function of the Semiconductor1D
"""
DistFun(S::Semiconductor1D,E) = DistFun(S.points,E)
"""
    C(S::Semiconductor1D) Signed doping concentration of the semiconductor
"""
C(S::Semiconductor1D) = C.(S.points)
"""
    φfₑ(S::Semiconductor1D) Electron Quasi-Fermi level
"""
φfₑ(S::Semiconductor1D) = φfₑ.(S.points)
"""
    φfₚ(S::Semiconductor1D) holes Quasi-Fermi level
"""
φfₚ(S::Semiconductor1D) = φfₚ.(S.points)

"""
    translate!(S::Semiconductor1D,pt::Number) Translates the semiconductor slab by pt 
"""
function translate!(S::Semiconductor1D,pt::Number)
    v = S.points
    for i in 1:length(v)
        S.points[i] = SemiconductorPoint1D(S[i],S[i].p+pt)
    end
end

function vjoint(L1D::Semiconductor1D,R1D::Semiconductor1D,dx=1e-9)
    if R₀(R1D) < Rₗ(L1D)
        translate!(R1D,Rₗ(L1D)+dx)
    end
    #Create Log Spaced Points on Transition
    lastL = pop!(L1D.points)
    firstR = popfirst!(R1D.points)
    midle = 0.5*(r(L1D[end])+r(R1D[1]))
    r1 = collect(10.0.^LinRange(log10(r(L1D[end])),log10(midle),11))
    popfirst!(r1)
    pop!(r1)
    #Left Side
    for r in r1
        push!(L1D,SemiconductorPoint1D(lastL.Data,r,ψ(lastL)))
    end


    r2 = collect(10.0.^LinRange(log10(midle),log10(r(R1D[1])),11))
    pop!(r2)
    popfirst!(r2)
    modelR = R1D[1]
    for r in reverse(r2)
        pushfirst!(R1D.points,SemiconductorPoint1D(firstR.Data,r,ψ(firstR)))
    end

    ΔE = (Ef(L1D[1])-Ef(R1D[1]))
    applyPsiProfile(R1D,ΔE.*ones(length(R1D)))
    #Δχ(R1D,-ΔE.*ones(length(R1D)))
    v = vcat(L1D.points,R1D.points)
    Semiconductor1D(v)
end

function vjointLin(L1D::Semiconductor1D,R1D::Semiconductor1D,dx=1e-9)
    if R₀(R1D) < Rₗ(L1D)
        translate!(R1D,Rₗ(L1D))
    end
    #Create Log Spaced Points on Transition
    lastL = pop!(L1D.points)
    #firstR = popfirst!(R1D.points)
    
    ΔE = (Ef(L1D[1])-Ef(R1D[1]))
    Δχ(R1D,ΔE)
    applyPsiProfile(R1D,-ΔE.*ones(length(R1D)))
    
    v = vcat(L1D.points,R1D.points)
    Semiconductor1D(v)
end


"""
    BandDiagram(S::Semiconductor1D) Plots the Semiconductor1D BandDiagram
"""
function BandDiagram(S::Semiconductor1D)
    BandDiagram(S,Rₗ(S))
end

function BandDiagram(S::Semiconductor1D,LM)
    y = 1e4*r(S)
	if(maximum(y) < LM)
        Sp = SemiconductorPoint(S[end].Data,LM,0.0)
		push!(S,Sp)
	end
    ψₛ = ψ(S[1])
    #plot(r(S),Ec(S)+χ(S),label=L"E_0")
    plot(r(S),Ec(S),label=L"E_C")
    plot!(r(S),Ev(S),label=L"E_V")
    plot!(r(S),Ef(S)-ψ(S),label=L"\varPhi+E_F")
    plot!(r(S),Ei(S),label=L"E_i")

    LM = L(S)
    return plot!(label=L"E_F",xlims=(0,LM),xlabel=L"r\, [\mu m]",legend=:best)#,yticks=nothing)

end

function BandDiagramJ(S::Semiconductor1D)
    BandDiagramJ(S,Rₗ(S))
end

function BandDiagramJ(S::Semiconductor1D,LM)
    #plot(r(S),Ec(S)+χ(S),label=L"E_0")
    plot(1e4*r(S),Ec(S),label=L"E_C")
    plot!(1e4*r(S),Ev(S),label=L"E_V")
    #plot!(1e4*r(S),Ei(S),label=L"\varPhi+E_F")
    plot!(1e4*r(S),-ψ(S)+Ei(S),label=L"E_i")

    #plot!(r(S),Ei(S),label=L"E_i")

    LM = L(S)*1e4
    return plot!(label=L"E_F",xlims=(0,LM),xlabel=L"r\, [\mu m]",legend=:best)#,yticks=nothing)

end

function CarrierInjectionPlot(S::Semiconductor1D,LM = 5.0)
    y = r(prob)
    LM = maximum(y)
    Ny = N(prob)
    Py = P(prob)
    niy = ni(prob)
    plot(y,Ny./niy,label=L"N/n_i")
    plot!(y,Py./niy,label=L"P/n_i",xlims=(-0.1,LM),yaxis=:log10,xlabel=L"r\, [\mu m]")
end


"""
    SingleSlabBVProblem(ψₛ,S::SemiconductorMaterial,f::DistributionFunc,T::Number) Solves by integrating the analyticall 1/ξ(x) the uniform semiconductor SingleSlabBVProblem and return a SemiconductorPoint1D
"""
function SingleSlabBVProblem(ψₛ,S::SemiconductorMaterial,f::DistributionFunc,T::Number)
    ## Precalculate the Material Properties
    Sd = SemiconductorData(T,f,S)
    if ψₛ == 0
        y = collect(LinRange(0,150,101))
        psi = 0*y
        return Semiconductor1D(Sd,y,psi)
    end
    psi1 = LinRange(ψₛ,0.5*ψₛ,51)
    psi2 = 10.0.^LinRange(log10(abs(ψₛ*0.5)),log10(abs(ψₛ*1e-3)),101)
    if ψₛ < 0
        psi2 = -1.0 .* psi2
    end
    psi = vcat(psi1,psi2)
    SdV = [SemiconductorData(Sd,p*Eunit) for p in psi]
    integrand_y(ψ) = 1/ξ(Sd,ψ)
    Y(ψ) = quadgk(integrand_y,ψ,ψₛ)[1]
    L = abs.(Y.(psi)) 
    return Semiconductor1D(return Semiconductor1D(Sd,L,psi))
end

function relaxationCostFunction(S::Semiconductor1D)
    Eest = -∇ρ.(Ref(S),r(S))
    error = sqrt(sum((ξ(S).-Eest).^2))
end

"""
    applyPsiProfile(S::Semiconductor1D,ψ::T) where T <: AbstractArray
    
    Apply an electrostatic potential to the semiconductor 1D
"""
function applyPsiProfile(S::Semiconductor1D,ψ::T) where T <: AbstractArray
    if length(S) != length(ψ)
        throw(DomainError(ψ,"ψ must have the same number of points as spatial points as S"))
    end
    for i in 1:length(ψ)
        S[i] = applyPsi(S[i],ψ[i])
    end
end

"""
    applyEfeProfile(S::Semiconductor1D,φfₑ::AbstractArray)
    
    Apply an electron quasi-fermi level to the semiconductor 1D
"""
function applyEfeProfile(S::Semiconductor1D,φfₑ::AbstractArray)
    if length(S) != length(φfₑ)
        throw(DomainError(ψ,"ψ must have the same number of points as spatial points as S"))
    end
    for i in 1:length(φfₑ)
        S[i] = applyEfe(S[i],φfₑ[i])
    end
end
"""
    applyEfpProfile(S::Semiconductor1D,φfₚ::AbstractArray)
    
    Apply an holes quasi-fermi level to the semiconductor 1D
"""
function applyEfpProfile(S::Semiconductor1D,φfₚ::AbstractArray)
    if length(S) != length(φfₚ)
        throw(DomainError(φfₚ,"ψ must have the same number of points as spatial points as S"))
    end
    for i in 1:length(φfₚ)
        S[i] = applyEfp(S[i],φfₚ[i])
    end
end

function ApplyPsiAtPoint(S::Semiconductor1D,ψ,x)
    v = r(S)
    if x in v
        i = findlast(v .== x)
        S[i] = applyPsi(S[i],ψ)
    else
        m = findlast( v .< x )
        M = m + 1
        δr = v[M]-v[m]
        δM = v[M]-x
        δm = x-v[m]
        S[M] = applyPsi(S[M],δM/δr *ψ)
        S[m] = applyPsi(S[m],δm/δr *ψ)
    end
end

