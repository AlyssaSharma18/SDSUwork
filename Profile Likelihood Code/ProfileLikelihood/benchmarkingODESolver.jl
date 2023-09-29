using ParameterizedFunctions, MATLABDiffEq, OrdinaryDiffEq,
      ODEInterfaceDiffEq, Plots, Sundials, SciPyDiffEq
using DiffEqDevTools
using LinearAlgebra, StaticArrays
using BenchmarkTools

f = @ode_def_bare SISVectorHostModel begin
  dI_h = beta_h * S_h * I_v - (mu_h + gamma) * I_h
  dI_v = beta_v * S_v * I_h - mu_v * I_v
  dS_h = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h
  dS_v = Pi_v - mu_v * S_v - beta_v * S_v * I_h
  dC_h = beta_h * S_h * I_v
end Pi_h Pi_v mu_h mu_v beta_h beta_v gamma
p = [0.004, 90, 0.00004, 0.09, 0.0001, 0.001, 0.09] 
tspan = (0.0, 50.0)
u0 = [1.0, 1.0, 1000.0, 5000.0, 1.0]
prob = ODEProblem(f,u0,tspan,p)

sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

# Finding overhead of MATLAB code 
mlAlg = MATLABDiffEq.ode45()
algstr = string(typeof(mlAlg).name.name)

# @benchmark solve($prob, $mlAlg)
# @benchmark MATLABDiffEq.eval_string("[t,u] = $(algstr)(diffeqf,tspan,u0,options);")

setups = [
          Dict(:alg=>DP5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>TsitPap8())
          Dict(:alg=>Vern7())
          Dict(:alg=>Vern9())
          Dict(:alg=>Feagin14())
          Dict(:alg=>MATLABDiffEq.ode23())
          Dict(:alg=>MATLABDiffEq.ode45())
          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>SciPyDiffEq.RK45())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.odeint())
  ]

labels = [
  "Julia: DP5"
  "Julia: Tsit5"
  "Julia: TsitPap8"
  "Julia: Vern7"
  "Julia: Vern9"
  "Julia: Feagin14"
  "MATLAB: ode23"
  "MATLAB: ode45"
  "MATLAB: ode113"
  "SciPy: RK45"
  "SciPy: LSODA"
  "SciPy: odeint"
  ]

abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet([prob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      appxsol=[test_sol],dense=false,
                      save_everystep=false,numruns=100,maxiters=10000000,
                      timeseries_errors=false,verbose=false)
plt = plot(wp,title="SIS Vector-Host Model",legend=:outertopleft,size=(800,500),
     color=permutedims([repeat([:LightGreen],6)...,
    repeat([:Orange],3)...,repeat([:Red],3)...]),
     xticks = 10.0 .^ (-12:1:5),
     yticks = 10.0 .^ (-6:0.5:5),
     bottom_margin=5Plots.mm,
     dpi = 800)
display(plt)
# savefig(plt, "benchmarkPlot.png")