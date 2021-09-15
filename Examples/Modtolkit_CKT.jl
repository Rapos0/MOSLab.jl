using ModelingToolkit, Plots, DifferentialEquations


rc_eqs = [
          connect(Vg.p, M1.g)
          connect(R1.n, M1.d,C1.p)
          connect(Vdd.p, R1.p)
          connect(M1.s,C1.n,Vdd.n,Vg.n,ground.g)


         ]

@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,
                          [R1, C1, M1,Vdd,Vg, ground])
sys = structural_simplify(rc_model)
u0 = [
      C1.v => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 0.01))
sol = solve(prob, Tsit5())
plot(sol)