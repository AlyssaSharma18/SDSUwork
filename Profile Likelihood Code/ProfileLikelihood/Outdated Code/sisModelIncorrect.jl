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
 tspan = (0., 40.)
 p0 = [0.0001, 0.001, 0.09, 0.00004, 0.09]
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
perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.01 * i), lower = 0), prob, times, solver_opts)

# Initial guess that was calculated from repeated trials of optimization
p0 = [0.0005518, 0.000261035, 0.505069, 2.71894e-16, 0.0718024]

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :de_rand_2_bin, 
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 1.,
    :maxTimeSecond => 1.,
    :traceMode => :silent
)

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p0, fitter_opts, noisyData, 1, prob, solver_opts, times)
println("The minimum loss is $loss.") 
println("The fitted parameters are $paramsFitted.")

# Find threshold 
threshold = ProfileLikelihood.findThreshold(.95, 5, loss)
thresholdP = ProfileLikelihood.findThreshold(.95, 1, loss)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :de_rand_2_bin, 
    :algSecond => :de_rand_2_bin_radiuslimited,
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 30.,
    :maxTimeSecond => 30.,
    :traceMode => :compact
)

# theta, sol = ProfileLikelihood.PL(0.0000025, 50, 1, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
# PLbeta_h = plot(theta, [sol, (x) -> threshold], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Threshold"], dpi = 400)
# scatter!([paramsFitted[1]], [loss], color = "green", labels = "Fitted Parameter")
# savefig(PLbeta_h, "PLbeta_h.png")

# CHANGE THIS TO GET A BETTER LOOKING GRAPH
theta3Improved, sol3Improved = ProfileLikelihood.PL(0.000003, 510, 4, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
PLmuh = plot(theta3Improved, [sol3Improved, (x) -> threshold], xlabel = L"\mu_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[4]], [loss], color = "green", labels = "Fitted Parameter")
savefig(PLmuh, "PLmu_h.png")

# theta1, sol1 = ProfileLikelihood.PL(0.000001, 100, 2, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
# PLbeta_v = plot(theta1, [sol1, (x) -> threshold], xlabel = L"\beta_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[2]], [loss], color = "green", labels = "Fitted Parameter")
# savefig(PLbeta_v, "PLbeta_v.png")

# theta2, sol2 = ProfileLikelihood.PL(0.00001, 100, 3, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
# PLgamma = plot(theta2, [sol2, (x) -> threshold], xlabel = L"\gamma", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[3]], [loss], color = "green", labels = "Fitted Parameter")
# savefig(PLgamma, "PLgamma.png")

# theta4, sol4 = ProfileLikelihood.PL(0.0002, 40, 5, paramsFitted, noisyData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL)
# PLmu_v = plot(theta4, [sol4, (x) -> threshold], xlabel = L"\mu_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Threshold"], right_margin=5mm, dpi = 400)
# scatter!([paramsFitted[5]], [loss], color = "green", labels = "Fitted Parameter")
# savefig(PLmu_v, "PLmu_v.png")

jldsave("sisModelData2.jld2"; theta3Improved, sol3Improved)

# sisModelData -> threshold, thresholdP, theta, sol, theta1, sol1, theta2, sol2, theta4, sol4 
# sisModelData2 -> theta3Improved, sol3Improved