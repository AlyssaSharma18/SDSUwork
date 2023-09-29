# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


# Define system of ODEs
# Constants 
const Pi_h = 0.004
const Pi_v = 90

# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, C_h, beta_h, beta_v, gamma, mu_h, mu_v) = (u[1], u[2], u[3], u[4], u[5], p[1], p[2], p[3], p[4], p[5])
        du[1] = beta_h * S_h * I_v - (mu_h + gamma) * I_h
        du[2] = beta_v * S_v * I_h - mu_v * I_v
        du[3] = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h
        du[4] = Pi_v - mu_v * S_v - beta_v * S_v * I_h
        du[5] = beta_h * S_h * I_v
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 40.0)
p0 = [0.0001, 0.001, 0.09, 0.00004, 0.09]
u0 = [1.0, 1.0, 1000.0, 5000.0, 1.0]

p1 = [0.0011, 0.011, 0.19, 0.00014, 0.19]

prob = ODEProblem(sis!, u0, tspan, p0);
prob1 = remake(prob, p = p1)

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-4,
    :abstol => 1e-6,
)

# times
times = LinRange{Float64}(0.0, 40.0, 50)

# Generate data for testing 
perfectData1, noisyData1 = ProfileLikelihood.generateData(5, 366, i -> truncated(Normal(i, 0.01 * i), lower=-eps(Float64)), prob, times, solver_opts; prevalenceStatus = true)

# Testing Likelihood function 
p1 = [0.0011, 0.011, 0.19, 0.00014, 0.19]
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.01)

# res1 = ProfileLikelihood.likelihood(p1, perfectData, 1, prob, solver_opts, times, obj)
# res2 = ProfileLikelihood.likelihood(p1, perfectData, 2, prob, solver_opts, times, obj)
# res3 = ProfileLikelihood.likelihood(p1, perfectData, 3, prob, solver_opts, times, obj)
# res4 = ProfileLikelihood.likelihood(p1, perfectData, 4, prob, solver_opts, times, obj)
# res5 = ProfileLikelihood.likelihood(p1, perfectData, 5, prob, solver_opts, times, obj)
# res = res1 + res2 

ProfileLikelihood.likelihood(p1, [perfectData1], [], prob, solver_opts, times, [obj]; prevalenceObserved = [5])

perfectData2, noisyData2 = ProfileLikelihood.generateData(5, 364, i -> truncated(Poisson(i), lower=-eps(Float64)), prob, times, solver_opts; prevalenceStatus = true)
println(perfectData2)

# Plotting stuff 
probCur = remake(prob, p=p0)
sol = solve(prob, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol],saveat=times)
# plt = plot(times, sol[5,:], legend=:topleft, labels = "Fitted Parameters", dpi = 400)
plt = plot(times, perfectData2, labels = "Perfect Data")
plot!(times, noisyData2, labels = "Noisy Data")
display(plt)