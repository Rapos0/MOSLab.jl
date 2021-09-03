abstract type BandGapModel end


"""
Experimental Data for the Silicon Band Gap on eV and Kelvin
"""
TexpBG = [3.902699812884201,19.24084469393206,75.93156909917133,88.55564465829099,110.23790430366208,168.89779916243427,192.92167869553595,248.05845139445782,289.0171968279426,360.29047491758,411.88274080014264]
BgexpD = [1.1557025624729165,1.1555143933922498,1.1530489952245262,1.1520632962588473,1.1492423449784364,1.1405983597348612,1.1351549126101401,1.1235500821966806,1.113266923299148,1.093717856420627,1.0793094811010167].+0.013587437527083468




struct VarshniModel <: BandGapModel
    p::Vector{Number}
    f::Function
end

struct GengModel <: BandGapModel 
    p::Vector{Number}
    f::Function
end

"""
# Arguments 
    T:: Number Latice Temperature [K]
    f:: GengModel
Calculate the Thermal Capacity at constant pressure for the GengBandGap Model
"""
function Cp(T::Number,m::GengModel)
    Θ = m.p[5]
    if(T <= 0)
        return 4*pi^4/15
    else
        G(x) = x^4*exp(x)/(exp(x)-1)^2
        return quadgk(G,0,Θ/T)[1]
    end
end

"""
# Arguments 
    p:: [Band Gap at T0[eV], T0 [K], Band Gap at T1[eV], T1[K]]
Calculation of Geng Band Gap Model for low doping Semicondurs
"""
function GengModel(p=[BgexpD[1],TexpBG[1],BgexpD[end],TexpBG[end],480])
    return GengModel(p,x-> x)
end


"""
# Arguments 
    p:: [α [eV/K²], β [K] , Eg(0) [eV]]
Calculation of Band Gap Temperature Dependence Based on Varshin Model Eg(T) = Eg(0) - αT²/(T+β)
"""
VarshniModel(;α=7.021e-4*Eunit/Tunit/Tunit,β=1108*Tunit,Eg0 =1.16929*Eunit) = VarshniModel([Eg0,α,β],(T,p)-> p[1] - p[2]*T^2/(T+p[3]))
#VarshniModel(p = [1.16929,4.9e-4,650]) = VarshniModel(p,(T,p)-> p[1] - p[2]*T^2/(T+p[3]))
Eg(T::Number,m::BandGapModel) = m.f(T,m.p)

function IntGen(T::Number,m::GengModel)
    Θ = m.p[5]
    T0 = m.p[2]
    if T == 0
        return 0.0
    else
        X = collect(LinRange(1.0,T,450))
        Y = (X/Θ).^3 .* Cp.(X,Ref(m))
        return NumericalIntegration.integrate(X,Y)
    end
end

function Eg(T::Number,m::GengModel)
    T = ustrip(T)
    Eg0 = m.p[1]
    T0 = m.p[2]
    Eg1 = m.p[3]
    T1 = m.p[4]
    N = IntGen(T,m)
    D = IntGen(T1,m)
    Egn = (Eg0 - (Eg0-Eg1)*N/D)*Eunit
    return Egn
end

