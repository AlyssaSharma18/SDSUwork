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

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.01 * i), lower=0), prob, times, solver_opts)

# Initial guess that was calculated from repeated trials of optimization
p1 = [0.00010329841041789106, 0.0009732057799278893, 0.09114563530043511, 4.897370054495707e-16, 0.09116629763383406]

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :de_rand_2_bin,
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 60.0,
    :maxTimeSecond => 60.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.relativeError(data, sol, 0.01)

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, noisyData, 1, prob, solver_opts, times, obj)
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")

# Find threshold 
threshold = ProfileLikelihood.findThresholdOptimization(0.95, 5, loss)
thresholdP = ProfileLikelihood.findThresholdOptimization(0.95, 1, loss)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :de_rand_2_bin,
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 120.0,
    :maxTimeSecond => 60.0,
    :traceMode => :compact
)

# Plot solution
# probCur = remake(prob, p=paramsFitted)
# sol = solve(probCur, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol], saveat=times)
# plot(sol)
# scatter!(times, noisyData)

# Constants to add back 
PLconst = ProfileLikelihood.likelihoodConst("relativeError"; noiseLevel = 0.01, times)

# Graph profile likelihood 
# theta, sol = ProfileLikelihood.PL(0.0000003, 100, 1, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, "relativeError", obj)
# sol = sol .+ PLconst 
# PLbeta_h = plot(theta, [sol, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], dpi = 400)
# scatter!([paramsFitted[1]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# savefig(PLbeta_h, "PLbeta_h.png")

# theta3, sol3 = ProfileLikelihood.PL(0.0000006, 600, 4, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, obj)
# sol3 = sol3 .+ PLconst 
# PLmuh = plot(theta3, [sol3, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\mu_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[4]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# savefig(PLmuh, "PLmu_h.png")

# theta1, sol1 = ProfileLikelihood.PL(0.0000012, 210, 2, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, obj)
# sol1 = sol1 .+ PLconst
# PLbeta_v = plot(theta1, [sol1, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[2]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# savefig(PLbeta_v, "PLbeta_v.png")

theta2, sol2 = ProfileLikelihood.PL(0.001, 15, 3, paramsFitted, noisyData, 1, threshold + 30, loss, prob, solver_opts, times, fitter_opts_PL, obj)
sol2 = sol2 .+ PLconst
PLgamma = plot(theta2, [sol2, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\gamma", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[3]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
savefig(PLgamma, "PLgamma.png")

# theta4, sol4 = ProfileLikelihood.PL(0.00005, 20, 5, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, obj)
# sol4 = sol4 .+ PLconst
# PLmu_v = plot(theta4, [sol4, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\mu_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[5]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
# savefig(PLmu_v, "PLmu_v.png")

# jldsave("sisCorrectData.jld2"; loss, paramsFitted, threshold, thresholdP, PLconst, theta, sol, theta1, sol1, theta2, sol2, theta3, sol3, theta4, sol4)
# jldsave("sisCorrectData2.jld2"; theta4, sol4)