# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, Random, Distributions, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, StaticArrays


# Define system of ODEs
# Constants 
const Pi_h = 0.004
const Pi_v = 90
const mu_h = 0.00004
const mu_v = 0.09

# Define SIS model
function sis_static(u, p, t)
    dI_h = p[1] * u[3] * u[2] - (mu_h + p[3]) * u[1]
    dI_v = p[2] * u[4] * u[1] - mu_v * u[2]
    dS_h = Pi_h - mu_h * u[3] - p[1] * u[3] * u[2] + p[3] * u[1]
    dS_v = Pi_v - mu_v * u[4] - p[2] * u[4] * u[1]
    dC_h = p[1] * u[3] * u[2]
    SA[dI_h, dI_v, dS_h, dS_v, dC_h]
end

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
p0 = [0.0001, 0.001, 0.09] 
u0 = SA[1.0, 1.0, 1000.0, 5000.0, 1.0]

u1 = [1.0, 1.0, 1000.0, 5000.0, 1.0]

prob = ODEProblem(sis_static, u0, tspan, p0);
prob1 = ODEProblem(sis!, u1, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-5,
    :abstol => 1e-10,
)

# times
times = LinRange{Float64}(0.0, 30.0, 31)
println("The amount of data is ", length(times))

solve(prob, Tsit5(), reltol = solver_opts[:reltol], abstol = solver_opts[:abstol])

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob1, times, solver_opts; incidenceStatus = true)

# Plot solution 
plt = plot(times, perfectData, dpi = 400, legend = :topleft, labels = "Perfect Incidence Data", xlabel = L"t")
plot!(times, noisyData, labels = "Noisy Incidence Data")
display(plt)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 15.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.3)

# True loss value 
trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [], prob1, solver_opts, times, [obj]; incidenceObserved = [5])
println("The loss with the true parameters is ", trueLoss)

# Initial guess (after repeated trials of optimization)
p1 = [4.909723455681667e-5, 0.0035588417145903826, 0.5925358671958921]

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, [noisyData], [], prob1, solver_opts, times, [obj]; incidenceObserved = [5])
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")


