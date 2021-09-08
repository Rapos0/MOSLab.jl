

struct CircuitComponent
    name::String
    component::Component
    pins::Dict{String,Int}
end

CircuitComponent(comp::TransistorModel) = CircuitComponent(name,comp,["d"=>1,"s"=>2,"g"=>3])

mutable struct Circuit
    components::AbstractVector{CircuitComponent}
    connection_matrix::AbstractMatrix{Bool}
    nodes::Dict{String,Int}
    conduction_matrix
end

add_node!(ckt::Circuit,node::Pair{String,Int}) = push!(ckt.nodes,node) 

Circuit() = Circuit([],zeros(Bool,1,1),Dict("gnd!"=> 0),[])

function add_component(ckt::Circuit,comp::TransistorModel,connections::Dict{String,Int})
    #Add Nodes to Circuit if they don't exists already
    for node in connections
        if node[2] ∉ ckt.nodes
            add_node!(ckt,node)
        end
    end

end
