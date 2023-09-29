# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, StatsPlots


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
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.3)

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

# Plot the four fitted parameters 
p5 = [5.2843886620615984e-5, 0.004130829427391291, 0.8698907140968367]
p15 = [5.1833216998212974e-5, 0.003958359788117503, 0.7857912956362569]
p25 = [5.0105966340021236e-5, 0.003697114411185875, 0.6589579952017401]
p35 = [4.8059906471622965e-5, 0.0034201038732419938, 0.5273912205564562]
p30 = [4.909724477651933e-5, 0.0035588322896510363, 0.592532802113042]
probCur5 = remake(prob, p=p5)
probCur15 = remake(prob, p=p15)
probCur25 = remake(prob, p=p25)
probCur35 = remake(prob, p=p35)
probCur30 = remake(prob, p=p30)
sol5 = ProfileLikelihood.generateIncidenceData(5, probCur5, times, solver_opts)
sol15 = ProfileLikelihood.generateIncidenceData(5, probCur15, times, solver_opts)
sol25 = ProfileLikelihood.generateIncidenceData(5, probCur25, times, solver_opts)
sol35 = ProfileLikelihood.generateIncidenceData(5, probCur35, times, solver_opts)
sol30 = ProfileLikelihood.generateIncidenceData(5, probCur30, times, solver_opts)
pltNoiseLevel = plot(times, [sol5, sol15, sol25, sol35, sol30], legend=:topright, labels = ["5% Noise level" "15% Noise Level" "25% Noise Level" "35% Noise Level" "30% Noise Level"], xlabel = L"t", dpi = 400)
scatter!(times, noisyData, labels = "Noise Incidence Data")
display(pltNoiseLevel)

relativeError5 = (abs.(sol5) - noisyData) ./ abs.(sol5)
relativeError15 = (abs.(sol15) - noisyData) ./ abs.(sol15)
relativeError25 = (abs.(sol25) - noisyData) ./ abs.(sol25)
relativeError35 = (abs.(sol35) - noisyData) ./ abs.(sol35)
relativeError30 = (abs.(sol30) - noisyData) ./ abs.(sol30)

deleteat!(relativeError5, 1)
deleteat!(relativeError15, 1)
deleteat!(relativeError25, 1)
deleteat!(relativeError35, 1)
deleteat!(relativeError30, 1)

d5 = fit(Normal{Float64}, relativeError5)
d15 = fit(Normal{Float64}, relativeError15)
d25 = fit(Normal{Float64}, relativeError25)
d35 = fit(Normal{Float64}, relativeError35)
d30 = fit(Normal{Float64}, relativeError30)

println(d5)
println(d15)
println(d25)
println(d35) 
println(d30)

# Going to choose a noise level of 30% 