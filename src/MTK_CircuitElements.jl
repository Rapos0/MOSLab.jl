@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, []; name=name)
end

function ModelingToolkit.connect(::Type{Pin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

function MKT_Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
end

function MKT_AbstractTransistor(;name)
    @named d = Pin()
    @named g = Pin()
    @named s = Pin()
    sts = @variables vds(t)=1.0 ids(t)=1.0 vgs(t)=1.0
    eqs = [
           vds ~ d.v - s.v
           vgs ~ g.v - s.v
           0 ~ d.i + s.i + g.i
           g.i ~ 0
           ids ~ d.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name), d, g,s)
end


function MKT_OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

function MKT_Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function MKT_Transistor(;name, W=1.0,L=1.0,m=QuadraticModel())
    @named tr = AbstractTransistor()
    @unpack vds, ids, vgs = tr
    ps = @parameters K=K Vt=Vt L=L
    eqs = [
           ids ~ K*(vgs-Vt)^2*(1+vds*L)
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), tr)
end

function MKT_Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function MKT_AcVoltage(;name, V = 1.0, Vac = 0.1, f=1e3)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           v ~ V + Vac*sin(2*π*f*t) 
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function MKT_ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end
