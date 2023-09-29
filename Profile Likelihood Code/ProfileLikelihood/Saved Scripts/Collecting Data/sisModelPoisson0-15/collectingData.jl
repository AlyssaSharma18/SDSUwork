# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics

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
    :searchRange => [(0.0, 0.2), (0.0, 0.5), (0.0, 2.0)],
    :maxTimeFirst => 180.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.poissonError(data, sol)

trueLossArr = Vector{Float64}()
bbLossArr = Vector{Float64}()
mhLossArr = Vector{Float64}()
bbParamsArr = Vector{Vector{Float64}}()
mhParamsArr = Vector{Vector{Float64}}()
close = 0
off = 0
# 1002:1130
for seed in 1:200
    perfectData, noisyData = ProfileLikelihood.generateData(5, seed, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

    trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
    println("The loss with the true parameters is ", trueLoss)
    push!(trueLossArr, trueLoss)
    p1 = [0.0001, 0.001, 0.09]
    loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
    println("The minimum loss is $loss.")
    println("The fitted parameters are $paramsFitted.")
    push!(bbLossArr, loss)
    push!(bbParamsArr, paramsFitted)
end
jldsave("collectionData.jld2"; bbLossArr, bbParamsArr, trueLossArr)

# TODO: Test code on actual likelihood and graph them to see what they look like, also idea is to provide an initla guess and hope that it does not bias it 