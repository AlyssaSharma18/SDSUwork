# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, StatsPlots, MATLAB


# Define system of ODEs
# Constants 
const Pi_h = 0.004
const Pi_v = 90
const mu_h = 0.00004
const mu_v = 0.09

# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, C_h, beta_h, beta_v, gamma) = (u[1], u[2], u[3], u[4], u[5], p[1], p[2], p[3])
        du[1] = beta_h * S_h * I_v - (mu_h + gamma) * I_h
        du[2] = beta_v * S_v * I_h - mu_v * I_v
        du[3] = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h
        du[4] = Pi_v - mu_v * S_v - beta_v * S_v * I_h
        du[5] = beta_h * S_h * I_v
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 50.0)
p0 = [0.0001, 0.001, 0.09] 
u0 = [1.0, 1.0, 1000.0, 5000.0, 1.0]

prob = ODEProblem(sis!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-5,
    :abstol => 1e-10,
)

# times
times = LinRange{Float64}(0.0, 30.0, 31)
println("The amount of data is ", length(times))

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

# Plot solution 
plt = plot(times, perfectData, dpi = 400, legend = :topleft, labels = "Perfect Incidence Data", xlabel = L"t")
plot!(times, noisyData, labels = "Noisy Incidence Data")
display(plt)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :generating_set_search, # adaptive_de_rand_1_bin_radiuslimited
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 200.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.constVarianceError(data, sol, 1.0)

# True loss value 
# trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
# println("The loss with the true parameters is ", trueLoss)

# Initial guess (after repeated trials of optimization)
# use the true parameter as a guess 

# Find optimal parameters 
# loss, paramsFitted = ProfileLikelihood.estimateParams(p0, fitter_opts, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
# println("The minimum loss is $loss.")
# println("The fitted parameters are $paramsFitted.")

# Plot fitted parameters 
# probCur = remake(prob, p=paramsFitted)
# sol1 = ProfileLikelihood.generateIncidenceData(5, probCur, times, solver_opts)
# plt1 = plot(times, sol1, legend=:topleft, labels = "Fitted Parameters", xlabel = L"t", dpi = 400)
# scatter!(times, noisyData, labels = "Noisy Incidence Data")
# display(plt1)

# Plot the one fitted parameters 
p1 = [8.699983460903588e-5, 0.0011597338843241933, 0.10880269944252362]
probCur5 = remake(prob, p=p1)
ol1 = ProfileLikelihood.generateIncidenceData(5, probCur5, times, solver_opts)
pltConstVar = plot(times, [sol1], legend=:topright, labels = ["Constant Variance"], xlabel = L"t", dpi = 400)
scatter!(times, noisyData, labels = "Noisy Incidence Data")
display(pltConstVar)

constVar1 = (abs.(sol1 - noisyData))

deleteat!(constVar1, 1)

d1 = fit(Normal{Float64}, constVar1)

mat"""
    [$d2, $d3] = mle($(constVar1),'pdf',@(x,sigma)normpdf(x,0,sigma),'start',1.5)
"""
println(d2)
println(d1)

# choose a sigma of 5.7