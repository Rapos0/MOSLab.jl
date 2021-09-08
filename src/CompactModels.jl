abstract type Component end

abstract type TransistorModel <: Component end

#= struct QuadraticModel <: TransistorModel
    parameters::Dict{Symbol,Any}
    material::SemiconductorMaterial
end =#

Id(Vg,Vd,Vs,Vb,m::TransistorModel) = 0.0

Id(Vg,Vd,Vs,m::TransistorModel) = Id(Vg,Vd,Vs,Vs,m)

gm(Vg,Vd,Vs,Vb,m::TransistorModel) = gradient(x-> Id(x, Vd, Vs, Vb,m),Vg)[1]

gm(Vg,Vd,Vs,m::TransistorModel) = gradient(x-> Id(x, Vd, Vs,m),Vg)[1]

gmId(Vg,Vd,Vs,Vb,m::TransistorModel) = gradient(x-> log(Id(x, Vd, Vs, Vb,m)),Vg)[1]

gmId(Vg,Vd,Vs,m::TransistorModel) = gradient(x-> log(Id(x, Vd, Vs,m)),Vg)[1]

gds(Vg,Vd,Vs,Vb,m::TransistorModel) = gradient(x-> Id(Vg, x, Vs, Vb,m),Vd)[1]

gds(Vg,Vd,Vs,m::TransistorModel) = gradient(x-> Id(Vg, x, Vs,m),Vd)[1]

gdsId(Vg,Vd,Vs,Vb,m::TransistorModel) = gradient(x-> log(Id(Vg, x, Vs, Vb,m)),Vd)[1]

gdsId(Vg,Vd,Vs,m::TransistorModel) = gradient(x-> log(Id(Vg, x, Vs,m)),Vd)[1]
