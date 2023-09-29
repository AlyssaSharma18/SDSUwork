# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DifferentialEquations, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Optimization, OptimizationPolyalgorithms, OptimizationOptimJL, SciMLSensitivity, Zygote, Plots

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

# Initial condition
u0 = [1.0, 1.0, 1000.0, 5000.0]

# Simulation interval and intermediary points
tspan = (0.0, 40.0)
tsteps = LinRange{Float64}(0.0, 40.0, 30)

# LV equation parameter. p = [α, β, δ, γ]
p = [0.0001, 0.001, 0.09, 0.00004, 0.09]

# Setup the ODE problem, then solve
prob = ODEProblem(sis!, u0, tspan, p)
sol = solve(prob, Tsit5())

solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-5,
    :abstol => 1e-10,
)

perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.01 * i), lower=0), prob, tsteps, solver_opts)

function loss(p)
  sol = solve(prob, Rosenbrock23(), p=p, saveat = tsteps)
  loss = sum(abs2, sol.-1)
  return loss, sol
end

# the likelihood function is not right and so, I need to fix it to make the code work 
function likelihood(paramsCur, data, prob)
    # solve odes
    sol = solve(
        prob,
        p=paramsCur,
        Rosenbrock23(),
        saveat=tsteps,
    )
    # loss
    # return 10000.0 * sum(((sol - data)./sol).^2) + 2 * sum(log.(sol))
    return sum((sol[1,:] - data) .^ 2), sol
    # return sum(sol) - sum(data.*log.(sol))
end


adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p)->likelihood(x, noisyData, prob), adtype)
optprob = Optimization.OptimizationProblem(optf, p)

result_ode = Optimization.solve(optprob, PolyOpt(), maxiters = 10000)