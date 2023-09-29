# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


# Define system of ODEs
# Define SIS model
function sis!(du, u, p, t)
    let (I_h, I_v, S_h, S_v, beta_h, beta_v, gamma, mu_h, mu_v, Pi_h, Pi_v) = (u[1], u[2], u[3], u[4], p[1], p[2], p[3], p[4], p[5], p[6], p[7])
        du[1] = beta_h * S_h * I_v - (mu_h + gamma) * I_h
        du[2] = beta_v * S_v * I_h - mu_v * I_v
        du[3] = Pi_h - mu_h * S_h - beta_h * S_h * I_v + gamma * I_h  
        du[4] = Pi_v - mu_v * S_v - beta_v * S_v * I_h
    end
end

 # Timespan, true parameters and initial conditions for simulating data
 tspan = (0., 40.)
 p0 = [0.0001, 0.001, 0.09, 0.00004, 0.09, 0.004, 90]
 u0 = [1., 1., 1000., 5000.]
 
 prob = ODEProblem(sis!, u0, tspan, p0);

 # solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-7,
    :abstol => 1e-14
)

# times
times = LinRange{Float64}(0., 40., 50)

# Generate data 
perfectData1, noisyData1 = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.1 * i), lower = 0), prob, times, solver_opts)
perfectData2, noisyData2 = ProfileLikelihood.generateData(2, 366, i -> truncated(Normal(i, 0.1 * i), lower = 0), prob, times, solver_opts)

noisyData = [noisyData1, noisyData2]

# Initial guess 
p1 = [0.00010315746323976762, 0.001019359727757813, 0.16367079985075036, 0.0015980817783814575, 0.09079870002009513, 13.782690452330739, 91.80466959088398]

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :de_rand_2_bin, 
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 1.0), (0.0, 2.0), (0.0, 2.0), (0.0, 1.0), (0.0, 2.0), (0.0, 50.0), (0.0, 200.0)],
    :maxTimeFirst => 10.,
    :maxTimeSecond => 5.,
    :traceMode => :compact
)

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, noisyData, [1 2], prob, solver_opts, times)
println("The minimum loss is $loss.") 
println("The fitted parameters are $paramsFitted.")

# Find threshold 
threshold = ProfileLikelihood.findThreshold(.95, 7, loss)
thresholdP = ProfileLikelihood.findThreshold(.95, 1, loss)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :de_rand_2_bin, 
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 1.0), (0.0, 2.0), (0.0, 2.0), (0.0, 1.0), (0.0, 2.0), (0.0, 50.0), (0.0, 200.0)],
    :maxTimeFirst => 30.,
    :maxTimeSecond => 15.,
    :traceMode => :compact
)

# theta1, sol1 = ProfileLikelihood.PL(0.0000001, 20, 1, paramsFitted, noisyData, [1 2], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
# PLbeta_h = plot(theta1, [sol1, (x) -> threshold, (x) -> thresholdP], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simulataneous Threshold" "Pointwise Threshold"], dpi = 400)
# scatter!([paramsFitted[1]], [loss], color = "green", labels = "Fitted Parameter")
# savefig(PLbeta_h, "PLbeta_h.png")

theta2, sol2 = ProfileLikelihood.PL(0.0000001, 20, 1, paramsFitted, noisyData, [1 2], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
PLbeta_v = plot(theta2, [sol2, (x) -> threshold, (x) -> thresholdP], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simulataneous Threshold" "Pointwise Threshold"], dpi = 400)
scatter!([paramsFitted[1]], [loss], color = "green", labels = "Fitted Parameter")
savefig(PLbeta_h, "PLbeta_h.png")