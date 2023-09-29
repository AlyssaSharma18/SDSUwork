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
times = LinRange{Float64}(0.0, 15.0, 16)
println("The amount of data is ", length(times))

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Normal(i, 0.15 * i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

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
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.35)

# True loss value 
# trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
# println("The loss with the true parameters is ", trueLoss)

# Initial guess (after repeated trials of optimization)
# will be the same as the true parameter 
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

# Plot the three fitted parameters 
p5 =  [0.00012793595194816588, 0.0010142873032018542, 0.3450634509772331]
p15 =  [0.0001254831975325942, 0.001009010429494083, 0.3260368641093373]
p25 = [0.00012110586501657963, 0.0009994870740920305, 0.29239207509439435]
p35 = [0.00011559290103614331, 0.0009868179438635707, 0.25025950922630724]
probCur5 = remake(prob, p=p5)
probCur15 = remake(prob, p=p15)
probCur25 = remake(prob, p=p25)
probCur35 = remake(prob, p=p35)
sol5 = ProfileLikelihood.generateIncidenceData(5, probCur5, times, solver_opts)
sol15 = ProfileLikelihood.generateIncidenceData(5, probCur15, times, solver_opts)
sol25 = ProfileLikelihood.generateIncidenceData(5, probCur25, times, solver_opts)
sol35 = ProfileLikelihood.generateIncidenceData(5, probCur35, times, solver_opts)
pltNoiseLevel = plot(times, [sol5, sol15, sol25, sol35], legend=:topright, labels = ["5% Noise level" "15% Noise Level" "25% Noise Level" "35% Noise Level"], xlabel = L"t", dpi = 400)
scatter!(times, noisyData, labels = "Noise Incidence Data")
display(pltNoiseLevel)

relativeError5 = (abs.(sol5) - noisyData) ./ abs.(sol5)
relativeError15 = (abs.(sol15) - noisyData) ./ abs.(sol15)
relativeError25 = (abs.(sol25) - noisyData) ./ abs.(sol25)
relativeError35 = (abs.(sol35) - noisyData) ./ abs.(sol35)

deleteat!(relativeError5, 1)
deleteat!(relativeError15, 1)
deleteat!(relativeError25, 1)
deleteat!(relativeError35, 1)

d5 = fit(Normal{Float64}, relativeError5)
d15 = fit(Normal{Float64}, relativeError15)
d25 = fit(Normal{Float64}, relativeError25)
d35 = fit(Normal{Float64}, relativeError35)

println(d5)
println(d15)
println(d25)
println(d35)

# Normal{Float64}(μ=0.00856830527182559, σ=0.10017469764509825)
# Normal{Float64}(μ=-0.010861329612780461, σ=0.10233046884203628)
# Normal{Float64}(μ=-0.04766792635109191, σ=0.10643973765769603)
# Normal{Float64}(μ=-0.0986726974975877, σ=0.11218748111422282)