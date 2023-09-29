# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, NLopt, FiniteDifferences


# Define system of ODEs
# Constants 
const Pi_h = 0.004
const Pi_v = 90

# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, beta_h, beta_v, gamma, mu_h, mu_v) = (u[1], u[2], u[3], u[4], p[1], p[2], p[3], p[4], p[5])
        du[1] = beta_h * S_h * I_v - (mu_h + gamma) * I_h
        du[2] = beta_v * S_v * I_h - mu_v * I_v
        du[3] = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h
        du[4] = Pi_v - mu_v * S_v - beta_v * S_v * I_h
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 40.0)
p0 = [0.0001, 0.001, 0.09, 0.00004, 0.09]
u0 = [1.0, 1.0, 1000.0, 5000.0]

prob = ODEProblem(sis!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-4,
    :abstol => 1e-6,
)

# times
times = LinRange{Float64}(0.0, 40.0, 50)

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.01 * i), lower=0), prob, times, solver_opts)

# Initial guess that was calculated from repeated trials of optimization
p1 = [0.00010329841041789106, 0.0009732057799278893, 0.09114563530043511, 4.897370054495707e-16, 0.09116629763383406]

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 7.0,
    :maxTimeSecond => 60.0,
    :traceMode => :compact,
    :algNLOpt => :LD_TNEWTON,
    :lowerBoundsNLOpt => [0., 0., 0., 0., 0.],
    :upperBoundsNLOpt => [2., 2., 2., 2., 2.]
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.01)
stepSize = 0.0002
g = (p0) -> ProfileLikelihood.likelihood(p0, noisyData, 1, prob, solver_opts, times, obj; paramIndex = 2, paramEval = 0.0009732057799278893 + stepSize)

h = (p0) -> ProfileLikelihood.likelihood(p0, noisyData, 1, prob, solver_opts, times, obj)

function q(p, gradient)
   if length(gradient) > 0
    gradient[:] = grad(central_fdm(4, 1, max_range=1e-4), g, p)[1]
  end
  return g(p)
end

function q1(p, gradient)
    if length(gradient) > 0
     gradient[:] = grad(central_fdm(4, 1, max_range=1e-4), h, p)[1]
   end
   return h(p)
 end

# p2 = [0.00010329841041789106, 0.0009732057799278893, 0.09114563530043511, 4.897370054495707e-16, 0.09116629763383406]
# println(q(p2, [1., 1., 1., 1., 1.]))
# println(q1(p2, [1., 1., 1., 1., 1.]))
p2 = [0.000103298, 0.0911456, 4.89737e-16, 0.0911663]
# println(q(p2, [1., 1., 1., 1., 1.]))

# res = q([0.0009732057799278893, 0.09114563530043511, 4.897370054495707e-16, 0.09116629763383406], [1., 1., 1., 1.])

println("LD_TNEWTON")
opt5 = Opt(:LD_TNEWTON, length(p2))
opt5.min_objective = (p, grad) -> q(p, grad)
opt5.lower_bounds = [0., 0., 0., 0.]
opt5.upper_bounds = [2., 2., 2., 2.]
opt5.ftol_rel = 1e-7
opt5.ftol_abs = 1e-14
(minf,minx,ret) = NLopt.optimize(opt5, p2)

resFirst = bboptimize(g, minx; Method=fitter_opts[:algFirst], SearchRange=fitter_opts[:searchRange], MaxTime=fitter_opts[:maxTimeFirst], TraceMode=fitter_opts[:traceMode], PopulationSize = 1000)

# Find optimal parameters 
# loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, noisyData, 1, prob, solver_opts, times, obj)
# println("The minimum loss is $loss.")
# println("The fitted parameters are $paramsFitted.")