# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, Random, Distributions, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


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
p0 = [0.0001, 0.001, 0.09] 
u0 = [1.0, 1.0, 1000.0, 5000.0, 1.0]

prob = ODEProblem(sis!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-5,  
    :abstol => 1e-10, 
)

# times
times = LinRange{Float64}(0.0, 30.0, 31)
println("The amount of data is ", length(times))

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

# Plot solution 
plt = plot(times, perfectData, dpi = 400, legend = :topleft, labels = "Perfect Incidence Data", xlabel = L"t")
plot!(times, noisyData, labels = "Noisy Incidence Data")
display(plt)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :generating_set_search, #adaptive_de_rand_1_bin_radiuslimited
    :searchRange => [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)],
    :maxTimeFirst => 20.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.constVarianceError(data, sol, 5.7)

# True loss value 
trueLoss = ProfileLikelihood.likelihood(p0, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
println("The loss with the true parameters is ", trueLoss)

# Initial guess (after repeated trials of optimization)
p1 = 

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, [noisyData], [], prob, solver_opts, times, [obj]; incidenceObserved = [5])
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")

# Plot fitted parameters 
probCur = remake(prob, p=paramsFitted)
# sol1 = solve(prob, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol],saveat=times)
sol1 = ProfileLikelihood.generateIncidenceData(5, probCur, times, solver_opts)
plt1 = plot(times, sol1, legend=:topleft, labels = "Fitted Parameters", xlabel = L"t", dpi = 400)
scatter!(times, noisyData, labels = "Noisy Incidence Data")
display(plt1)

# Find threshold 
threshold = ProfileLikelihood.findThresholdOptimization(0.95, 3, loss)
thresholdP = ProfileLikelihood.findThresholdOptimization(0.95, 1, loss)
println(threshold)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(eps(Float64), 0.5), (eps(Float64), 0.5), (eps(Float64), 40.)], 
    :maxTimeFirst => 60.0,
    :traceMode => :compact,
    :algNLOpt => :LD_TNEWTON,
    :lowerBoundsNLOpt => [eps(Float64), eps(Float64), eps(Float64)],
    :upperBoundsNLOpt => [1000., 1000., 1000.],
    :boundsMH => [eps(Float64) eps(Float64) eps(Float64); 0.5 0.5 40.]
)

# Constants to add back 
PLconst = ProfileLikelihood.likelihoodConst("constVarianceError"; times = times, sigma = 5.7)
println("The profile likelihood constant is ", PLconst)

[8.699928280353557e-5, 0.0011597447431600193, 0.10880420170900695]# beta_h 
theta1, sol1 = ProfileLikelihood.PL(8.4e-7, 55, 1, paramsFitted, [noisyData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], status = "globalBB")
sol1 = sol1 .+ PLconst 
PLbeta_h = plot(theta1, [sol1, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_h", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[1]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")

display(PLbeta_h)
# savefig(PLbeta_h, "PLbeta_h.png")

# beta_v 
theta2, sol2 = ProfileLikelihood.PL(1.27e-4, 95, 2, paramsFitted, [noisyData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], status = "globalBB")
sol2 = sol2 .+ PLconst 
PLbeta_v = plot(theta2, [sol2, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\beta_v", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[2]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLbeta_v)
# savefig(PLbeta_v, "PLbeta_v.png")

# gamma
theta3, sol3 = ProfileLikelihood.PL(0.04, 100, 3, paramsFitted, [noisyData], [], threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, [obj]; incidenceObserved = [5], status = "globalBB")
sol3 = sol3 .+ PLconst 
PLgamma = plot(theta3, [sol3, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\gamma", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[3]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")

display(PLgamma)
# savefig(PLgamma, "PLgamma.png")

# jldsave("sisIncidencePoissonMisfitConstantError0_30.jld2"; loss, paramsFitted, threshold, thresholdP, PLconst, theta1, sol1, theta2, sol2, theta3, sol3)