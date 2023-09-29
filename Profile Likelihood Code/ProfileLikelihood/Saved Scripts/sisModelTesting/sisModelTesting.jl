# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


# Define system of ODEs
# Constants 
const beta_v = 1/1000
const mu_h = 4/100000
const Pi_h = 1/250
const gamma = 9/100
const mu_v = 9/100

# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, Pi_v, beta_h) = (u[1], u[2], u[3], u[4], p[1], p[2])
        du[1] = beta_h * S_h * I_v - (mu_h + gamma) * I_h
        du[2] = beta_v * S_v * I_h - mu_v * I_v
        du[3] = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h
        du[4] = Pi_v - mu_v * S_v - beta_v * S_v * I_h
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 40.0)
p0 = [90.0, 0.0001]
u0 = [1.0, 1.0, 1000.0, 5000.0]

prob = ODEProblem(sis!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-7,
    :abstol => 1e-14,
)

# times
times = LinRange{Float64}(0.0, 40.0, 50)

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 1), lower=0), prob, times, solver_opts)

# Plot solution 
sol = solve(prob, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol],saveat=times)
plt = plot(times, sol[1,:], dpi = 400, legend = :topleft, labels = L"I_h", xlabel = L"t")
plot!(times, perfectData, labels = L"I_h")
display(plt)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :adaptive_de_rand_1_bin_radiuslimited,
    :searchRange => [(0.0, 1000.0), (0.0, 2.0)],
    :maxTimeFirst => 500.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.constVarianceError(data, sol, 1.0)

# True loss value 
trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [1], prob, solver_opts, times, [obj])
println("The loss with the true parameters is ", trueLoss)

# Initial guess 
p1 = [2.0, 2.0]

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, [noisyData], [1], prob, solver_opts, times, [obj])
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")

# Find threshold 
# threshold = ProfileLikelihood.findThresholdOptimization(0.95, 2, loss)
# thresholdP = ProfileLikelihood.findThresholdOptimization(0.95, 1, loss)
# println(threshold)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(0.0, 500.0), (0.0, 2.0)], 
    :maxTimeFirst => 5.0,
    :traceMode => :compact,
    :algNLOpt => :LD_TNEWTON,
    :lowerBoundsNLOpt => [0., 0.],
    :upperBoundsNLOpt => [500., 2.]
)

# Constants to add back 
# PLconst = ProfileLikelihood.likelihoodConst("constVarianceError"; times = times, sigma = 1.0)
# println(PLconst)

# gamma
# theta1, sol1 = ProfileLikelihood.PL(0.01, 7, 1, paramsFitted, [noisyData], [1], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj])
# println(sol1)
# sol1 = sol1 .+ PLconst 
# PLPi_v = plot(theta1, [sol1, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\Pi_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[1]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# display(PLPi_v)
# savefig(PLPi_v, "PLPi_v.png")

# beta_h
# theta2, sol2 = ProfileLikelihood.PL(0.00000001, 7, 2, paramsFitted, [noisyData], [1], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj])
# sol2 = sol2 .+ PLconst 
# PLbeta_h = plot(theta2, [sol2, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[2]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# display(PLbeta_h)
# savefig(PLbeta_h, "PLbeta_h.png")

# jldsave("sisTesting.jld2"; loss, paramsFitted, threshold, thresholdP, PLconst, theta1, sol1, theta2, sol2, theta3, sol3)




