struct Circuit{T}
    netlist::Netlist
    conduction_matrix::AbstractMatrix{T}
    currents_matrix#::Union{AbstractVector{T},SparseVector{T,Int},SparseVector{Number,Int}}
    parameters
 end

 function Circuit(net::Netlist,s::Num)
    N = nnodes(net)
    M = spzeros(N,N)
    I = spzeros(N)
    Ps = [s]
    for component in net.components
        m,i,p = SymStamp(component,N,s)
        M += m
        I += i
        push!(Ps,p...)
    end

    Gm = M[2:end,2:end]
    Im = I[2:end]
    #println(Ps)
    Circuit(net,Gm,Im,Ps)
 end

 function Circuit(net::Netlist)
    N = nnodes(net)
    M = spzeros(N,N)
    I = spzeros(N)
    Ps = []
    for component in net.components
        m,i,p = SymStamp(component,N)
        M += m
        I += i
        push!(Ps,p...)
    end

    Gm = M[2:end,2:end]
    Im = I[2:end]
    #println(Ps)
    Circuit(net,Gm,Im,Ps)
end

function CircuitFunction(ckt::Circuit)
    N = nnodes(ckt.netlist)-1 #exclude ground
    Vvector = @variables V[1:N]
    Vvector  = Symbolics.scalarize(Vvector)
    M = copy(ckt.conduction_matrix)
    I = copy(ckt.currents_matrix)
    S = size(M)
    for comp in ckt.netlist.components
        for l in 1:S[1]
            I[l] = substitute(I[l],symSubs(comp,Vvector[1]))
            for c in 1:S[2]
                M[l,c] = substitute(M[l,c],symSubs(comp,Vvector[1]))
            end
        end
    end
    return Matrix(M),Vector(I),Vector(Vvector[1])
end

function circuitFuntion(ckt::Circuit)
    Mm,Im,v = CircuitFunction(ckt)
    expr = Mm*v-Im
    cktfunc = eval(build_function(expr,v)[1])
    return  cktfunc#nlsolve(cktfunc,rand(nnodes(ckt.netlist-1))
    
end

struct CircuitProperties{T}
    ckt::Circuit
    Vi::Int
    Vo::Int
    Z::AbstractMatrix{T}
    Y::AbstractMatrix{T}
    T::AbstractMatrix{T}
    H::AbstractMatrix{T}
    VoltageGain::T
    TransconductanceGain::T
    Rin::T
    Rout::T
end


function CircuitProperties(ckt::Circuit,Vi::Int,Vo::Int)
    Z = inv(Matrix(ckt.conduction_matrix))
    Z = Z[[Vi,Vo],[Vi,Vo]]
    Y = Z
    try
        Y = inv(Z)
    catch
        Y = Matrix(ckt.conduction_matrix)[[Vi,Vo],[Vi,Vo]]
    end
    # T Matrix
    y11 = Y[1,1]
    y12 = Y[1,2]
    y21 = Y[2,1]
    y22 = Y[2,2]
    A = -y22/y21
    B = (-1/y21)
    C = ((y12*y21-y11*y22)/y21)
    D = (-y11/y21)  
    T = [A B;C D]
    #H Matrix

    A = 1/y11
    B = (-y12/y11)
    C = y21/y11
    D = ((y11*y22-y12*y21)/y11) 
    H = [A B;C D]

    VoltageGain = simplify(1/T[1,1])
    TransconductanceGain = simplify(1/T[1,2])
    Rin = simplify(1/Y[1,1])
    Rout = simplify(Z[2,2])

    CircuitProperties(ckt,Vi,Vo,Z,Y,T,H,VoltageGain,TransconductanceGain,Rin,Rout)
end