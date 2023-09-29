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
tspan = (0.0, 50.0)
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
times = LinRange{Float64}(0.0, 45.0, 46)
println("The amount of data is ", length(times))

# Generate data 
perfectData, noisyData = ProfileLikelihood.generateData(5, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

# Plot fitted parameters 
pRelativeError = [4.909723455681667e-5, 0.0035588417145903826, 0.5925358671958921]
pConstVar = [8.699928280353557e-5, 0.0011597447431600193, 0.10880420170900695]
pTrue = [0.0001, 0.001, 0.09] 

probRelativeError = remake(prob, p=pRelativeError)
probConstVar = remake(prob, p=pConstVar)
probTrue = remake(prob, p=pTrue)

solRelativeError = ProfileLikelihood.generateIncidenceData(5, probRelativeError, times, solver_opts)
solConstVar = ProfileLikelihood.generateIncidenceData(5, probConstVar, times, solver_opts)
solTrue = ProfileLikelihood.generateIncidenceData(5, probTrue, times, solver_opts)

plt1 = plot(times, solRelativeError, legend=:topright, labels = "Relative Error", xlabel = L"t", ylabel = "Incidence",dpi = 400, guidefontfamily = "Computer Modern")
plot!(times, solConstVar, legend=:topright, labels = "Constant Variance Error", xlabel = L"t", dpi = 400)
plot!(times, solTrue, legend=:topright, labels = "True Parameters", xlabel = L"t", dpi = 400)

scatter!(times, noisyData, labels = "Noisy Incidence Data", color = 5)
display(plt1)
# savefig(plt1, "varyingRelativeErrorConstVarGraph")

R_0(beta_v1, beta_h1, mu_h1, gamma1, mu_v1) = sqrt((beta_h1 * beta_v1 * 1000 * 5000)/((mu_h1 + gamma1) * mu_v1))

g(beta_v1, betah_1, gamma1) = R_0(beta_v1, betah_1, 0.00004, gamma1,0.09)

println("True: ", g(p0[2], p0[1], p0[3]))
println("Relative Error: ", g(pRelativeError[2], pRelativeError[1], pRelativeError[3]))
println("Constant Variance Error: ", g(pConstVar[2], pConstVar[1], pConstVar[3]))

trueR_0 = g(p0[2], p0[1], p0[3])
R_016 =  g(pRelativeError[2], pRelativeError[1], pRelativeError[3])
R_031 = g(pConstVar[2], pConstVar[1], pConstVar[3])

println((R_016 - trueR_0)/(trueR_0))
println((R_031 - trueR_0)/(trueR_0))