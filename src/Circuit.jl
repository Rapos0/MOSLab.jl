

struct CircuitComponent
    component::Component
    pins::Dict{String,Int}
end

CircuitComponent(comp::TransistorModel) = CircuitComponent(comp,["d"=>1,"s"=>2,"g"=>3])

mutable struct Circuit
    components::AbstractVector{CircuitComponent}
    connection_matrix::AbstractMatrix{Bool}
    nodes::Dict{String,Int}
    conduction_matrix
end

add_node!(ckt::Circuit,node::Pair{String,Int}) = push!(ckt.nodes,node) 

Circuit() = Circuit([],[[]],[],[])

# function add_component(ckt::Circuit,comp::TransistorModel,connections::Dict{String,Int})
#     for node in connections
#         if node[2] ∉ ckt.nodes
#             add_node!(ckt,node)
#     end
# end
