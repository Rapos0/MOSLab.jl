module MOSLab

using ExportAll
    using Unitful
    using QuadGK
    using LineSearches
    using Optim
    using NumericalIntegration
    using NLsolve
    using LinearAlgebra
    using DiffEqOperators
    #using Plots
    using Roots
    using LaTeXStrings
    using DifferentialEquations
    using DomainSets
    using DomainSets:×
    using NonlinearSolve
    using Zygote
    using ModelingToolkit
    using DiffRules
    using VoronoiFVM
    using Cubature
    using LinearAlgebra
    using ExtendableGrids
    using Symbolics
    using SparseArrays
    #using Latexify

    include("SemiconductorConstants.jl")
    include("BandGapModels.jl")
    include("Distributions.jl")

    DiffRules.@define_diffrule MOSLab.trapFermi(x) = :(DtrapFermi($x))
    @register trapFermi(x)
    
    
    include("Dopping.jl")
    include("Semiconductor.jl")
    include("SemiconductorData.jl")
    struct SiO2 <: AbstractMaterial
        χ::Number
        Ec::Number
        Ef::Number
        Ev::Number
        ϵr::Number
    end
    
    SiO2() = SiO2(0.95,4.5*Eunit,0*Eunit,-4.5*Eunit,3.9)
    Ec(M::AbstractMaterial) = M.Ec
    Ef(M::AbstractMaterial) = M.Ef
    Ev(M::AbstractMaterial) = M.Ev
    
    struct  Metal <: AbstractMaterial 
        φₘ::Number
    end
    
    φₘ(m::Metal,T) = m.φₘ
    
    struct NPoly <: AbstractMaterial
        Nb::Number
    end
    function φₘ(m::NPoly,T)
        SD = SemiconductorData(T,BoltzmanDist(),NSilicon(m.Nb))
        χ(SD)-Ec(SD)
    end
    
    Alluminium() = Metal(4.1*Eunit)
    
    struct MOSStructure
        metal::AbstractMaterial
        oxide::AbstractMaterial
        semiconductor::SemiconductorData
        tox::Number
        tsi::Number
    end
    include("PSP.jl")
    include("SemiconductorPoint.jl")
    include("1DSemiconductor.jl")
    include("MOS.jl")
    include("PN.jl")
    include("CompactModels.jl")
    include("PaoSah.jl")
    include("Brews.jl")
    include("UICM.jl")
    include("BSIM6.jl")
    include("EKV.jl")
    include("QuadradicModel.jl")
    include("MOSFET.jl")
    include("Solvers.jl")

    # Symbolics.@register gm(vg,vd,vs,m::TransistorModel)

    # Symbolics.@register Id(vg,vd,vs,m::TransistorModel)

    # Symbolics.@register IdNR(vg,vd,vs,m::TransistorModel)

    # Symbolics.@register gds(vg,vd,vs,m::TransistorModel)


    include("MTK_CircuitElements.jl")
    include("Netlist.jl")
    include("Stamps.jl")
    include("Circuit.jl")


    @exportAll()

end # module
