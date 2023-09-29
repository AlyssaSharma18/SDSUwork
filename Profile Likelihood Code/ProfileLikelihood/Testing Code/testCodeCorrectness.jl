# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2


# Define system of ODEs

# Define SIS model
function test!(du, u, p, t)
    let (a, b, c, d) = (p[1], p[2], p[3], p[4])
        du[1] = (a + b + c + d)
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 40.0)
p0 = [1., 2., 3., 4.]
u0 = [50.0]

prob = ODEProblem(test!, u0, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => Tsit5(),
    :reltol => 1e-4,
    :abstol => 1e-6,
)

# times
times = LinRange{Float64}(0.0, 40.0, 40)

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(1, 366, i -> truncated(Normal(i, 0.01 * i), lower=0), prob, times, solver_opts)

# Settings for optimization 
fitter_opts = Dict(
    :algFirst => :de_rand_2_bin_radiuslimited,
    :searchRange => [(-50.0, 50.0), (-50.0, 50.0), (-50.0, 50.0), (-50.0, 50.0)],
    :maxTimeFirst => 10.0,
    :traceMode => :compact
)

# Objective function 
obj = (data, sol) -> ProfileLikelihood.constVarianceError(data, sol, 1)

# Parameter guess 
p1 = [2., 3., 4., 5.]

# Find optimal parameters 
loss, paramsFitted = ProfileLikelihood.estimateParams(p1, fitter_opts, perfectData, 1, prob, solver_opts, times, obj)
println("The minimum loss is $loss.")
println("The fitted parameters are $paramsFitted.")

# Find threshold 
threshold = ProfileLikelihood.findThresholdOptimization(0.95, 4, loss)
thresholdP = ProfileLikelihood.findThresholdOptimization(0.95, 1, loss)

# Find profile likelihood 
fitter_opts_PL = Dict(
    :algFirst => :generating_set_search,
    :searchRange => [(-50.0, 50.0), (-50.0, 50.0), (-50.0, 50.0), (-50.0, 50.0)],
    :maxTimeFirst => 10.0,
    :traceMode => :compact,
    :algNLOpt => :LD_TNEWTON,
    :lowerBoundsNLOpt => [-50., -50., -50., -50.],
    :upperBoundsNLOpt => [50., 50., 50., 50.]
)

# Plot solution
probCur = remake(prob, p=paramsFitted)
sol = solve(probCur, solver_opts[:alg], reltol=solver_opts[:reltol], abstol=solver_opts[:abstol], saveat=times)
plt = plot(sol)
scatter!(times, perfectData)
display(plt)

# Constants to add back 
PLconst = ProfileLikelihood.likelihoodConst("relativeError"; noiseLevel = 0.01, times)

# Plots
theta, sol = ProfileLikelihood.PL(0.1, 50, 1, paramsFitted, perfectData, 1, threshold + 3, loss, prob, solver_opts, times, fitter_opts_PL, obj)
sol = sol .+ PLconst 
PLa = plot(theta, [sol, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"a", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], dpi = 400)
scatter!([paramsFitted[1]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLa)