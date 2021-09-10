
ConductanceStamp(G,pa::Int,pb::Int,N::Int) = (sparse([pa,pa,pb,pb],[pa,pb,pa,pb],[G,-G,-G,G],N,N),spzeros(N))
TransConductanceStamp(G,oa::Int,ob::Int,ia::Int,ib::Int,N::Int) = (sparse([oa,oa,ob,ob],[ia,ib,ia,ib],[G,-G,-G,G],N,N),spzeros(N))
VoltageSourceStamp(V,p::Int,m::Int,C::Int,N::Int) = (sparse([p,m,C,C],[C,C,p,m],[1,-1,-1,1],N,N),sparse([C],[1],[-V],N,1))
CurrentSourceStamp(I,p::Int,m::Int,N::Int) = (spzeros(N,N),sparse([p,m],[1,1],[-I,I],N,1))


function SymStamp(T::CircuitComponent{W},N::Int) where W <: TransistorModel
    name = T.name
    gmSym = Symbol("gm_{$name}")
    gdsSym = (Symbol("gds_{$name}"))
    IdSym = (Symbol("Id_{$name}"))
    syms = @variables  $gmSym $gdsSym $IdSym
    d = T.connections["d"]+1
    g = T.connections["g"]+1
    s = T.connections["s"]+1
    M1,I1 = TransConductanceStamp(syms[1],d,s,g,s,N)
    M2,I2 = ConductanceStamp(syms[2],d,s,N)
    #M3,I3 = ConductanceStamp(1/syms[1],g,s,N)
    M3,I3 = CurrentSourceStamp(syms[3],d,s,N)
    M4,I4 = CurrentSourceStamp(syms[1],d,s,N)
    addSymParameter.(Ref(T),syms)
    return M1+M2+M3,I1+I2+I3,syms
end

function SymStamp(T::CircuitComponent{W},N::Int,p::Num) where W <: TransistorModel
    name = T.name
    gmSym = Symbol("gm_{$name}")
    gdsSym = Symbol("gds_{$name}")
    CgsSym = Symbol("Cgs_{$name}")
    CgdSym = Symbol("Cgd_{$name}")
    #CdsSym = Symbol("Cds_{$name}")
    syms = @variables  $gmSym $gdsSym $CgsSym $CgdSym 
    d = T.connections["d"]+1
    g = T.connections["g"]+1
    s = T.connections["s"]+1
    M1,I1 = TransConductanceStamp(syms[1],d,s,g,s,N)
    M2,I2 = ConductanceStamp(syms[2],d,s,N)
    M3,I3 = ConductanceStamp(1/syms[1],g,s,N)
    M4,I4 = ConductanceStamp(p*syms[3],g,s,N)
    M5,I5 = ConductanceStamp(p*syms[4],g,d,N)
    addSymParameter.(Ref(T),syms)
    #M6,I6 = ConductanceStamp(p*syms[5],d,s,N)
    return M1+M2+M3+M4+M5,I1+I2+I3+I4+I5,syms
end


function SymStamp(T::CircuitComponent{Resistor},N::Int,s=0)
    name = T.name
    G = Symbol("$name")
    syms = @variables $G
    p = T.connections["p"]+1
    n = T.connections["n"]+1
    M2,I2 = ConductanceStamp(1/syms[1],p,n,N)
    addSymParameter.(Ref(T),syms)
    return M2,I2,syms
end


function SymStamp(T::CircuitComponent{VoltageSource},N::Int,s=0)
    name = T.name
    G = Symbol("$name")
    syms = @variables $G
    p = T.connections["p"]+1
    n = T.connections["n"]+1
    j = T.connections["j"]+1
    addSymParameter.(Ref(T),syms)
    M2,I2 = VoltageSourceStamp(syms[1],p,n,j,N)
    return M2,I2,syms
end


function SymStamp(T::CircuitComponent{CurrentSource},N::Int,s=0)
    name = T.name
    G = Symbol("$name")
    syms = @variables $G
    p = T.connections["p"]+1
    n = T.connections["n"]+1
    addSymParameter.(Ref(T),syms)
    M2,I2 = CurrentSourceStamp(syms[1],p,n,N)
    return M2,I2,syms
end

function SymStamp(T::CircuitComponent{TransConductance},N::Int,s=0)
    name = T.name
    G = Symbol("$name")
    syms = @variables $G
    ip = T.connections["ip"]+1
    in = T.connections["im"]+1
    op = T.connections["op"]+1
    on = T.connections["om"]+1
    addSymParameter.(Ref(T),syms)
    M2,I2 = TransConductanceStamp(syms[1],op,on,ip,in,N)
    return M2,I2,syms
end