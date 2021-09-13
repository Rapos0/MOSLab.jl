using MOSLab
using Plots
using Revise
using ModelingToolkit
using NonlinearSolve
using DiffEqFlux
using Measurements
using Distributions
#= eqs = Mf(V)*V-If(V) .~ 0.0
@named ns = NonlinearSystem(eqs, V, [])
prob = NonlinearProblem(ns,guess,[])
sol = solve(prob) =#
@parameters Ra Rb

R1 = CircuitComponent("R1",Resistor(Ra))
R2 = CircuitComponent("R2",Resistor(Rb))
V1 = CircuitComponent("Vi",VoltageSource(5.0))
Rcnet = Netlist()
addComponent(Rcnet,R1,Dict(["p"=>2,"n"=>1]))
addComponent(Rcnet,R2,Dict(["p"=>1,"n"=>0]))
addComponent(Rcnet,V1,Dict(["p"=>2,"n"=>0]))

ckt = Circuit(Rcnet)



M,IM,V = CircuitFunction(ckt)
p = [Ra,Rb]
Mf = eval(build_function(M,V,p)[1])
If = eval(build_function(IM,V,p)[1])
eqs = 0 .~ Mf(V,p)*V-If(V,p)
@named ns = NonlinearSystem(eqs,V,p)
pv = [1e3,1e3]
prob = NonlinearProblem(ns,rand(3),[Rb=>pv[2] ± 0.1*pv[2],Ra=>pv[1] ± 0.1*pv[1]])
V = solve(prob)
pstd(V[1])
function loss(ns,p)
    prob = NonlinearProblem(ns,rand(3),[Rb=>p[2] ± 0.1*p[2],Ra=>p[1] ± 0.1*p[1]])
    V = solve(prob)
    return pstd(V[1])
end
result = DiffEqFlux.sciml_train(p->loss(ns,p),[1e3,1e3];maxiters=1000,lb=[50],ub=[1e6])
plot(10.0.^(2:0.2:6),p->loss(ns,[p,p]))
plot!(xaxis=:log10)