# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, StatsPlots, LatinHypercubeSampling

f = jldopen("testingDataCollection.jld2", "r")
bbLossArr = read(f, "bbLossArr")
trueLossArr = read(f, "trueLossArr")
mhLossArr = read(f, "mhLossArr")
bbParamsArr = read(f, "bbParamsArr")
mhParamsArr = read(f, "mhParamsArr")
closeOkay = read(f,"close")
close(f)

overfitArr = trueLossArr - bbLossArr
histogram(overfitArr, normalize = :probability)
plot!(Chisq(3))

indexArr = Vector{Int64}()
bbLossArrInterest = Vector{Float64}()
for (ind, val) in enumerate(bbLossArr)
    if abs(val - mhLossArr[ind]) > 0.001 
      push!(indexArr, ind)
      push!(bbLossArrInterest, val)
    end
end 

# println(indexArr)

indexArr = [1, 7, 18, 26, 36, 45, 47, 54, 62, 70, 80, 97, 100, 103, 115, 128]
# println("Off attempts: ", length(indexArr))
# println("Total attempts: ", length(bbLossArr))
# for ind in indexArr
#     println("Index:", ind)
#     println("BlackBoxOptim: ", bbLossArr[ind])
#     println("Metaheuristics: ", mhLossArr[ind])
#     println("BlackBoxOptim: ", bbParamsArr[ind])
#     println("Metaheuristics: ", mhParamsArr[ind])
# end 

indexArr = indexArr .+ 1001
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
tspan = (0.0, 40.0)
p0 = [0.0001, 0.001, 0.09] # TODO: Change these parameters 
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

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :adaptive_de_rand_1_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 220.0,
    :traceMode => :silent
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.poissonError(data, sol)

# num = 17
# println(indexArr[num])
# println(bbLossArrInterest[num])
closeCounter = 0
for (ind, seed) in enumerate(1002:1005)
    perfectData, noisyData = ProfileLikelihood.generateData(5, seed, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

    g(l) = ProfileLikelihood.likelihood(l, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5]) 

    bounds = [0.000001 0.000001 0.000001; 2. 2. 2.]

    options = Options(time_limit = 60.0, iterations = 1000000, f_calls_limit = 10000000000)
    information = Information()
    algor = DE(CR = 0.01, strategy = :rand1, options = options, information = information)
    # algor = ECA(η_max = 1.0, adaptive = true, resize_population = true, options = options, information = information)

    resOpt = optimize(g, bounds, algor)
    println(resOpt)
    fx = minimum(resOpt)
    println(fx)
    if abs(bbLossArr[ind] - fx) < 0.001
        global closeCounter += 1
    end
end

println("Current: ", closeCounter)
println("To beat: ", closeOkay)