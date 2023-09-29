# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, ColorSchemes

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

# Measurement Error with mean 0 and standard deviation 10
plt1 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> Normal(i, 10), prob, times, solver_opts)
    plot!(plt1, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt1)

# Measurement Error with mean 0 and standard deviation 0.01 y(t_i)
plt2 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> Normal(i, 0.01 * i), prob, times, solver_opts)
    plot!(plt2, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt2)

# Data with Poisson Distribution of Mean y(t_i)
plt3 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> Poisson(i), prob, times, solver_opts)
    plot!(plt3, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt3)

# Measurement Error with mean 0 and standard deviation 0.1 y(t_i)
plt4 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> Normal(i, 0.1 * i), prob, times, solver_opts)
    plot!(plt4, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt4)

# Measurement Error with mean 0 and standard deviation 50
plt5 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> Normal(i, 50), prob, times, solver_opts)
    plot!(plt5, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt5)

# Measurement Error with mean 0 and standard deviation 0.1/0.2/0.3 y(t_i) (truncated)
plt9 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 0.1 * i), lower = 0), prob, times, solver_opts)
    plot!(plt9, noisyDataPlot, labels = "", palette = palette(:linear_kgy_5_95_c69_n256, 100))

    _, noisyDataPlot2 = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 0.2 * i), lower = 0), prob, times, solver_opts)
    plot!(plt9, noisyDataPlot2, labels = "", palette = palette(:linear_kry_0_97_c73_n256, 100))

    _, noisyDataPlot1 = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 0.3 * i), lower = 0), prob, times, solver_opts)
    plot!(plt9, noisyDataPlot1, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt9)

# Measurement Error with mean 0 and standard deviation 50/100 (truncated)
plt13 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 10), lower = 0), prob, times, solver_opts)
    plot!(plt13, noisyDataPlot, labels = "", palette = palette(:linear_kgy_5_95_c69_n256, 100))
    _, noisyDataPlot1 = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 50), lower = 0), prob, times, solver_opts)
    plot!(plt13, noisyDataPlot1, labels = "", palette = palette(:linear_kry_0_97_c73_n256, 100))
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, 100), lower = 0), prob, times, solver_opts)
    plot!(plt13, noisyDataPlot, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt13)

# Measurement Error normally distributed with mean 0 and standard deviation with square root of mean and the poisson distrubution and both truncated 
plt14 = plot(xlabel = "Time", ylabel = L"I_h", dpi = 400)
for seed in 1:100
    _, noisyDataPlot = ProfileLikelihood.generateData(1, seed, i -> truncated(Normal(i, sqrt(i)), lower = 0), prob, times, solver_opts)
    plot!(plt14, noisyDataPlot, labels = "", palette = palette(:linear_kry_0_97_c73_n256, 100))
    _, noisyDataPlot1 = ProfileLikelihood.generateData(1, seed, i -> truncated(Poisson(i), lower = 0), prob, times, solver_opts)
    plot!(plt14, noisyDataPlot1, labels = "", palette = palette(:linear_blue_95_50_c20_n256, 100))
end 
display(plt14)
