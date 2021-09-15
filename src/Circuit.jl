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

function SymbolicMatrixes(ckt::Circuit)
    N = nnodes(ckt.netlist)-1
    Vvector = @variables V[1:N]
    Vvector  = Symbolics.scalarize(Vvector)
    M = copy(ckt.conduction_matrix)
    I = copy(ckt.currents_matrix)
    return Matrix(M),Vector(I),Vector(Vvector[1])
end

function SymbolicAuxEqs(ckt::Circuit,V)
    eq = []
    for comp in ckt.netlist.components
        dd = symSubs(comp,V)
        for d in dd
            push!(eq,first(d)~last(d))
        end
    end
    return eq
end


function CircuitFunction(ckt::Circuit)
    N = nnodes(ckt.netlist)-1 #exclude ground
    Vnames = map(x->x.name,ckt.netlist.nodes)[2:end]
    test = Symbolics.Symbol.(Vnames)
    Vvector = Vector(Symbolics.Variable.(test))
    M = copy(ckt.conduction_matrix)
    I = copy(ckt.currents_matrix)
    S = size(M)
    for comp in ckt.netlist.components
        for l in 1:S[1]
            I[l] = substitute(I[l],symSubs(comp,Vvector))
            for c in 1:S[2]
                M[l,c] = substitute(M[l,c],symSubs(comp,Vvector))
            end
        end
    end
    return Matrix(M),Vector(I),Vector(Vvector)
end


function dc_op(ckt::Circuit; maxiter=1000,abs_tol=1e-6,rel_tol=1e-2)
    Mm,Im,Vm = CircuitFunction(ckt)

    M_f = eval(build_function(Mm,Vm)[2])
    I_f = eval(build_function(Im,Vm)[2])

    v₀ = zeros(length(Vm))
    Mm = zeros(size(Mm))
    Im = zeros(size(Im))
    Base.invokelatest(M_f,Mm,v₀)
    Base.invokelatest(I_f,Im,v₀)
    v₁ = Mm\Im
    iter = 0
    abserror = 100
    relerror = 100
    while iter < maxiter && abserror > abs_tol && relerror > rel_tol
        v₀ = v₁
        Base.invokelatest(M_f,Mm,v₀)
        Base.invokelatest(I_f,Im,v₀)
        v₁ = Mm\Im
        iter += 1
        abserror = abs(maximum(v₁-v₀))
        relerror = abs(maximum((v₁-v₀)/v₁))
        if iter == maxiter
            @warn "Maximum iteration achieved"
        end
    end 
    return v₁   
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