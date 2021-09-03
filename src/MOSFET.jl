struct SemiconductorInterval{T}
    semiconductor::SemiconductorData
    S::T
end

function belongs(x,S::SemiconductorInterval )
    DomainSets.SVector(x...) in S.S
end

#∈(x,S::SemiconductorInterval) = belongs(x,S)

const SemiconductorDataInput = AbstractVector{SemiconductorInterval}

struct MOSFETInputDeck{T}
    Source::SemiconductorInterval
    Gate::MOSStructure
    Drain::SemiconductorInterval
    region
    DrainContact::T
    GateContact::T
    SourceContact::T
end
getSemi(x,S) = S[belongs.(Ref(x),S)][1].semiconductor

function MOSFETInputDeck(N_d,Nb,T,l,h,xd,ld,tox,type;lscal=1.0,ϕₘ=4.4,gate_ov=nothing)
    if type == :PMOS
        SD = SemiconductorData(T,BoltzmanDist(),PSilicon(N_d))
        SG = SemiconductorData(T,BoltzmanDist(),NSilicon(Nb))
    else 
        SD = SemiconductorData(T,BoltzmanDist(),NSilicon(N_d))
        SG = SemiconductorData(T,BoltzmanDist(),PSilicon(Nb))
    end 
    L = 2*ld+l
    if isnothing(lscal)
        lscal = L
    end
    if isnothing(gate_ov)
        gate_ov = 0.1*ld/lscal
    else
        gate_ov = gate_ov/lscal
    end
    SourceInt = (0..ld/lscal) × ((h-xd)/lscal..h/lscal)
    DrainInt = ((l+ld)/lscal..(L)/lscal) × ((h-xd)/lscal..h/lscal)
    GateInt = (0..L/lscal) × (0..h/lscal)
    GateMOS = MOSStructure(Metal(ϕₘ),SiO2(),SG,tox,h)
    Source = SemiconductorInterval(SD,SourceInt)
    Drain = SemiconductorInterval(SD,DrainInt)
    SourceContact = (0.0..0.8*ld/lscal)
    GateContact = (ld/lscal-gate_ov..(ld+l)/lscal+gate_ov)
    DrainContact = ((ld+l)/lscal+0.2*ld/lscal..L/lscal)
    MOSFETInputDeck(Source,GateMOS,Drain,GateInt,DrainContact,GateContact,SourceContact)
end

function C(x,y,fetIn::MOSFETInputDeck)
    Cv = [C(fetIn.Source.semiconductor),C(fetIn.Gate.semiconductor),C(fetIn.Drain.semiconductor)]
    return dot(Cv,belongVector(x,y,fetIn))
end

L(fetIN::MOSFETInputDeck) = fetIN.region.domains[1].right
H(fetIN::MOSFETInputDeck) = fetIN.region.domains[2].right

function belongVector(x,y,fetIn::MOSFETInputDeck)
    if DomainSets.SVector((x,y)) ∈ fetIn.region
        if  belongs((x,y),fetIn.Source)
            return [1,0,0]
        elseif belongs((x,y),fetIn.Drain)
            return [0,0,1]
        else
            return [0,1,0]
        end
    else
        println("($(x),$(y)) does not belong to the Grid")
        return [0.0,0.0,0.0]
    end
end

function Grid(fetIn::MOSFETInputDeck,nx,ny,lscal)
    l = L(fetIn)/lscal
    h = H(fetIn)/lscal
    x = LinRange(0.0,l,nx)
    y = LinRange(0.0,h,ny)
    dx = diff(x)[1]
    dy = diff(y)[1]
    X = VoronoiFVM.Grid(x,y)
    VoronoiFVM.bfacemask!(X, [(fetIn.SourceContact.left)/lscal,h-dy/2.0],[fetIn.SourceContact.right/lscal,h+dy/2.0],6)
    VoronoiFVM.bfacemask!(X, [fetIn.GateContact.left/lscal,h-dy/2.0],[fetIn.GateContact.right/lscal,h+dy/2.0],7)
    VoronoiFVM.bfacemask!(X, [fetIn.DrainContact.left/lscal,h-dy/2.0],[fetIn.DrainContact.right/lscal,h+dy/2.0],8)

    return X

end

function MOSFETSimulation(fetIn::MOSFETInputDeck,Vg,Vd,nx,ny;verbose=false,ξ₀=1e-10)
    Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(fetIn.Gate.semiconductor,L(fetIn))
    Rf(n,p) = 0.0
    Cf(x,y) = C(x,y,fetIn)
    Jn = q*1400*Nv(fetIn.Source.semiconductor)*Vscal/Jscal
    Jp = q*734*Nc(fetIn.Source.semiconductor)*Vscal/Jscal
    G = Grid(fetIn,nx,ny,lscal)

    #u,nn,pp = DDVoronoiMOSFET(G,Ec(fetIn.Source.semiconductor),Ev(fetIn.Source.semiconductor),Nc(fetIn.Source.semiconductor),Nv(fetIn.Source.semiconductor),Ec(fetIn.Source.semiconductor)-Ec(fetIn.Gate.semiconductor),Vg,0.0,lscal,Rscal,Cscal,Vscal,λ,Jn,Jp,Rf,Cf,fetIn;verbose=verbose,ξ₀=ξ₀)
    u,nn,pp = DDVoronoiMOSFET(G,Ec(fetIn.Source.semiconductor),Ev(fetIn.Source.semiconductor),Nc(fetIn.Source.semiconductor),Nv(fetIn.Source.semiconductor),Ec(fetIn.Source.semiconductor)-Ec(fetIn.Gate.semiconductor),Vg,Vd,lscal,Rscal,Cscal,Vscal,λ,Jn,Jp,Rf,Cf,fetIn;verbose=verbose,ξ₀=ξ₀)
    return G[Coordinates],u,nn,pp
end

function MOSFETSimulation(N_d,Nb,T,l,h,Vg,Vd,Nx,Ny)
    SD = SemiconductorData(T,BoltzmanDist(),NSilicon(N_d))
    SG = SemiconductorData(T,BoltzmanDist(),PSilicon(Nb))
    #= SourceInt = (0..l/3.0) × (h/2.0..h)
    Source = SemiconductorInverval(SD,SourceInt)
    DrainInt = (2.0*l/3.0..l) × (h/2.0..h) 
    Drain = SemiconductorInverval(SD,DrainInt)
    GateInt = setdiff((0..l)×(0..h),SourceInt)
    GateInt = setdiff(GateInt,DrainInt)
    Gate = SemiconductorInverval(SG,GateInt) =#

    #InputDeck = [Drain,Gate,Source]


    Vscal,lscal,Cscal,μscal,Dscal,Jscal,Rscal,tscal,δ,λ,Ld0,n_i = scalling(SG,l)
    Rf(n,p) = 0.0
    function Cf(x,y) 
        if x < l/3.0 && y > 0.5h
            return C(SD)
        elseif x > 2.0l/3.0 && y > 0.5h
            return C(SD)
        else
            return C(SG)
        end
    end
    Plots.contour(LinRange(0,l,101),LinRange(0,h,101),Cf,fill=true)

    Jn = q*1400*Nv(SD)*Vscal/Jscal
    Jp = q*734*Nc(SD)*Vscal/Jscal

    x = LinRange(0,l/lscal,Nx)
    dx = diff(x)[1]
    y = LinRange(0,h/lscal,Ny)
    dy = diff(y)[1]
    X = VoronoiFVM.Grid(x,y)
    0.2l/lscal-dx
    1.1*h/lscal
    VoronoiFVM.bfacemask!(X, [0.0,h/lscal-dy/2.0],[0.2l/lscal,h/lscal+dy/2.0],6)
    VoronoiFVM.bfacemask!(X, [0.25*l/lscal,h/lscal-dy/2.0],[0.75l/lscal,h/lscal+dy/2.0],7)
    VoronoiFVM.bfacemask!(X, [0.8*l/lscal,h/lscal-dy/2.0],[1.1l/lscal,h/lscal+dy/2.0],8)
    
    u,nn,pp = DDVoronoiMOSFET(X,Ec(SD),Ev(SD),Nc(SD),Nv(SD),Ec(SD)-Ec(SG),Vg,0.0,lscal,Rscal,Cscal,Vscal,Ld0,l,λ,21,Jn,Jp,21,Rf,Cf)

    u,nn,pp = DDVoronoiMOSFET(X,Ec(SD),Ev(SD),Nc(SD),Nv(SD),Ec(SD)-Ec(SG),Vg,Vd,lscal,Rscal,Cscal,Vscal,Ld0,l,λ,21,Jn,Jp,21,Rf,Cf;u0=u/Vscal,n0=nn,p0=pp)
    Cords = X[Coordinates]
    y1 = unique(Cords[2,:])[2]
    mask = Cords[2,:] .== y1
    u0 = u[mask]
    u0 .-= u0[1]
    return (Cords,Ec(SD).-u0,nn*Cscal,pp*Cscal)
end
