

struct CircuitComponent
    name::String
    component::Component
    pins::AbstractVector{String}
end

CircuitComponent(name::String,comp::TransistorModel) = CircuitComponent(name,comp,["d","g","s"])

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

add_node!(net::Netlist,node::Pair{String,Int}) = push!(net.nodes,node) 
add_node!(net::Netlist,i::Int) = push!(net.nodes,Node(i)) 

function addComponent(net::Netlist,comp::CircuitComponent,con::Dict{String,Int})
    push!(net.components,comp)
    for c in con
        println(first(c))
        if last(c) ∉ pos.(net.nodes)
            add_node!(net,last(c))
        end
    end
end




