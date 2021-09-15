mutable struct Resistor <: Component
    R
end

mutable struct VoltageSource <: Component
    V
end

mutable struct CurrentSource <: Component
    I
end

mutable struct TransConductance <: Component
    G
end

mutable struct Capacitor <: Component
    C
end

mutable struct CircuitComponent{T<: Component}
    name::String
    component::T
    connections::Dict{String,Union{Int,Nothing}}
    symParametes::AbstractVector
end

pins(circomp::CircuitComponent) = collect(keys(circomp.connections))
CircuitComponent(name::String,comp::TransistorModel) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("d"=>nothing,"g"=>nothing,"s"=>nothing),[])
CircuitComponent(name::String,comp::Resistor) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("p"=>nothing,"n"=>nothing),[])
CircuitComponent(name::String,comp::VoltageSource) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("p"=>nothing,"n"=>nothing,"j"=>nothing),[])
CircuitComponent(name::String,comp::CurrentSource) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("p"=>nothing,"n"=>nothing),[])
CircuitComponent(name::String,comp::TransConductance) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("op"=>nothing,"om"=>nothing,"ip"=>nothing,"im"=>nothing),[])
CircuitComponent(name::String,comp::Capacitor) = CircuitComponent(name,comp,Dict{String,Union{Int,Nothing}}("p"=>nothing,"n"=>nothing),[])

function symSubs(T::CircuitComponent{W},v) where W <: TransistorModel
    Vg = T.connections["g"] == 0 ? 0.0 : v[T.connections["g"]]
    Vd = T.connections["d"] == 0 ? 0.0 : v[T.connections["d"]]
    Vs = T.connections["s"] == 0 ? 0.0 : v[T.connections["s"]]
    d = Dict([T.symParametes[1]=> gm(Vg,Vd,Vs,T.component), T.symParametes[2] => gds(Vg,Vd,Vs,T.component),T.symParametes[3] => Id(Vg,Vd,Vs,T.component)-gm(Vg,Vd,Vs,T.component)*(Vg-Vs)-gds(Vg,Vd,Vs,T.component)*(Vd-Vs)gm(Vg,Vd,Vs,T.component)*(Vg-Vs)-gds(Vg,Vd,Vs,T.component)*(Vd-Vs)])
    #println(d)
    return d
end

function symSubs(T::CircuitComponent{Capacitor},v)
    return Dict([T.symParametes[1]=> T.component.C])
end

function symSubs(T::CircuitComponent{Resistor},v)
    return Dict([T.symParametes[1]=> T.component.R])
end

function symSubs(T::CircuitComponent{VoltageSource},v)
    return Dict([T.symParametes[1]=> T.component.V])
end

function symSubs(T::CircuitComponent{CurrentSource},v)
    return Dict([T.symParametes[1]=> T.component.I])
end

function symSubs(T::CircuitComponent{TransConductance},v)
    return Dict([T.symParametes[1]=> T.component.G])
end

addSymParameter(comp::CircuitComponent,p::Num) = push!(comp.symParametes,p)



struct Node
    name::String
    position::Int
end

Node(i::Int) = Node("V_n$(i)",i)
pos(n::Node) = n.position

mutable struct Netlist
    components::AbstractVector{CircuitComponent}
    nodes::AbstractVector{Node}
end

nnodes(net::Netlist) = length(net.nodes)

Netlist() = Netlist([],[Node("gnd!",0)])



add_connection(cc::CircuitComponent,pin::String,node::Int) = cc.connections[pin] = node 

add_node!(net::Netlist,node::Node) = push!(net.nodes,node) 
add_node!(net::Netlist,i::Int) = push!(net.nodes,Node(i)) 

function addComponent(net::Netlist,comp::CircuitComponent,con::Dict{String,Int})
    for (i,c) in enumerate(con)
        if last(c) ∉ pos.(net.nodes)
            add_node!(net,last(c))
        end
        add_connection(comp,first(c),last(c))
    end
    push!(net.components,comp)
end


function addComponent(net::Netlist,comp::CircuitComponent{VoltageSource},con::Dict{String,Int})
    for (i,c) in enumerate(con)
        if last(c) ∉ pos.(net.nodes)
            add_node!(net,last(c))
        end
        add_connection(comp,first(c),last(c))
    end
    push!(net.components,comp)
    Jnode = nnodes(net)
    add_node!(net,Node("J_$(comp.name)",Jnode))
    add_connection(comp,"j",Jnode)
end



