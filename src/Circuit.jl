

mutable struct CircuitComponent{T<: Component}
    name::String
    component::T
    connections::Dict{String,Union{UInt,Nothing}}
end

pins(circomp::CircuitComponent) = collect(keys(circomp.connections))
CircuitComponent(name::String,comp::TransistorModel) = CircuitComponent(name,comp,Dict{String,Union{UInt,Nothing}}("d"=>nothing,"g"=>nothing,"s"=>nothing))

Stamp(comp::CircuitComponent{TransistorModel}) = 


struct Node
    name::String
    position::Int
end

Node(i::Int) = Node("n$(i)",i)
pos(n::Node) = n.position

mutable struct Netlist
    components::AbstractVector{CircuitComponent}
    nodes::AbstractVector{Node}
end

nnodes(net::Netlist) = length(net.nodes)

Netlist() = Netlist([],[Node("gnd!",0)])


mutable struct Circuit
   netlist::Netlist
   conduction_matrix
end

add_connection(cc::CircuitComponent,pin::String,node::Int) = cc.connections[pin] = node 

add_node!(net::Netlist,node::Pair{String,Int}) = push!(net.nodes,node) 
add_node!(net::Netlist,i::Int) = push!(net.nodes,Node(i)) 

function addComponent(net::Netlist,comp::CircuitComponent,con::AbstractVector{Int})
    for (i,c) in enumerate(con)
        if last(c) ∉ pos.(net.nodes)
            add_node!(net,last(c))
        end
        add_connection(comp,pins(comp)[i],c)
    end
    push!(net.components,comp)
end




