# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


# Define system of ODEs
# Constants 
const Pi_h = 1.0/25000.0
const Pi_v = 1.0/10.0
const mu_h = 1.0/25000.0
const mu_v = 1.0/10.0

# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, C_h, gamma, beta_h, beta_v) = (u[1], u[2], u[3], u[4], u[5], p[1], p[2], p[3])
        du[1] = (beta_h * S_h * I_v) / (S_v + I_v) - (mu_h + gamma) * I_h
        du[2] = (beta_v * S_v * I_h) / (S_h + I_h) - (mu_v * I_v) 
        du[3] = Pi_h - (mu_h * S_h) - (beta_h * S_h * I_v) / (S_v + I_v) + gamma * I_h
        du[4] = Pi_v - (mu_v * S_v) - (beta_v * S_v * I_h) / (S_h + I_h) 
        du[5] = (beta_h * S_h * I_v) / (S_v + I_v)
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 216.0)
p0 = [0.1, 0.1, 0.2]
u0 = [0.01, 0.01, 0.99, 0.99, 0.01]

prob = ODEProblem(sis!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-7,
    :abstol => 1e-14,
)

# times
times = LinRange{Float64}(0.0, 216.0, 217)

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Normal(i, 0.01 * i), lower=-eps(Float64), upper = 1), prob, times, solver_opts; incidenceStatus = true)

# Plot solution 
sol = solve(prob, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol],saveat=times)
plt = plot(sol, dpi = 400, legend = :topleft, labels = [L"I_h" L"I_v" L"S_h" L"S_v" L"C_h"], xlabel = L"t")
plot!(times, perfectData, labels = "Prevalence Data")
display(plt)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :adaptive_de_rand_1_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 15.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.constVarianceError(data, sol, 0.001)

# True loss value 
trueLoss = ProfileLikelihood.likelihood(p0, [perfectData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
println("The loss with the true parameters is ", trueLoss)

# Initial guess 
p1 = [0.09999999999991685, 0.09999999999999966, 0.19999999999990375]

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, [perfectData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")

# Plot fitted parameters 
probCur = remake(prob, p=paramsFitted)
sol1 = ProfileLikelihood.generateIncidenceData(5, probCur, times, solver_opts)
plt1 = plot(times, sol1, legend=:topleft, labels = "Fitted Parameters", xlabel = L"t", ylabel = "Prevalence", dpi = 400)
scatter!(times, perfectData, labels = "Perfect Data")
display(plt1)

# Find threshold 
threshold = ProfileLikelihood.findThresholdOptimization(0.95, 3, loss)
thresholdP = ProfileLikelihood.findThresholdOptimization(0.95, 1, loss)
println(threshold)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(0.0, 5.0), (0.0, 5.0), (0.0, 5.0)], # for the profile likelihood of mu_h, the first parameter needed to be increased or it approaches the bounds 
    :maxTimeFirst => 25.0,
    :traceMode => :compact,
    :algNLOpt => :LD_TNEWTON,
    :lowerBoundsNLOpt => [0., 0., 0.],
    :upperBoundsNLOpt => [5., 5., 5.]
)

# Constants to add back 
PLconst = ProfileLikelihood.likelihoodConst("constVarianceError"; times = times, sigma = 0.001)
println(PLconst)

# gamma
theta1, sol1 = ProfileLikelihood.PL(5e-4, 2000, 1, paramsFitted, [perfectData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], order = 4, maxRange = 1e-6)
println(sol1)
sol1 = sol1 .+ PLconst 
PLgamma = plot(theta1, [sol1, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\gamma", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[1]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLgamma)
savefig(PLgamma, "PLgamma.png")

# beta_h
theta2, sol2 = ProfileLikelihood.PL(1e-3, 2000, 2, paramsFitted, [perfectData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], order = 5, maxRange = 1e-6)
sol2 = sol2 .+ PLconst 
PLbeta_h = plot(theta2, [sol2, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[2]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLbeta_h)
savefig(PLbeta_h, "PLbeta_h.png")

# beta_v
theta3, sol3 = ProfileLikelihood.PL(5e-4, 1000, 3, paramsFitted, [perfectData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], order = 5, maxRange = 1e-6)
sol3 = sol3 .+ PLconst 
PLmu_h = plot(theta3, [sol3, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\mu_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[3]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLmu_h)
savefig(PLmu_h, "PLmu_h.png")

jldsave("sisIncidenceData.jld2"; loss, paramsFitted, threshold, thresholdP, PLconst, theta1, sol1, theta2, sol2, theta3, sol3)




