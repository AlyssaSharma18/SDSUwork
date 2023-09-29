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
        du[6] = beta_v * S_v * I_h
    end
end

# Timespan, true parameters and initial conditions for simulating data
tspan = (0.0, 50.0)
p0 = [0.0001, 0.001, 0.09] 
u0 = [1.0, 1.0, 1000.0, 5000.0, 1.0, 1.0]

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

perfectDataHost, noisyDataHost = ProfileLikelihood.generateData(5, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)
perfectDataVector, noisyDataVector = ProfileLikelihood.generateData(6, 366, i -> truncated(Poisson(i), lower = -eps(Float64)), prob, times, solver_opts; incidenceStatus = true)

# Plot fitted parameters 
p15 = [9.592909229946155e-5, 0.001010520005939762, 0.0788359708967769]
p30 = [9.329466810800791e-5, 0.0010666198022741937, 0.09775869214627872]
pTrue = [0.0001, 0.001, 0.09] 

prob15 = remake(prob, p=p15)
prob30 = remake(prob, p=p30)
probTrue = remake(prob, p=pTrue)

sol15Host= ProfileLikelihood.generateIncidenceData(5, prob15, times, solver_opts)
sol30Host = ProfileLikelihood.generateIncidenceData(5, prob30, times, solver_opts)
solTrueHost = ProfileLikelihood.generateIncidenceData(5, probTrue, times, solver_opts)

sol15Vector = ProfileLikelihood.generateIncidenceData(6, prob15, times, solver_opts)
sol30Vector = ProfileLikelihood.generateIncidenceData(6, prob30, times, solver_opts)
solTrueVector = ProfileLikelihood.generateIncidenceData(6, probTrue, times, solver_opts)

plt1 = plot(times, sol15Host, legend=:topright, labels = "16 Data Points", xlabel = L"t", ylabel = "Incidence (Host)",dpi = 400, guidefontfamily = "Computer Modern")
plot!(times, sol30Host, legend=:topright, labels = "31 Data Points", xlabel = L"t", dpi = 400)
plot!(times, solTrueHost, legend=:topright, labels = "True Parameters", xlabel = L"t", dpi = 400)
scatter!(times, noisyDataHost, labels = "Noisy Incidence Host Data")
display(plt1)
# savefig(plt1, "hostIncidenceFittedTrue")

plt2 = plot(times, sol15Vector, legend=:topright, labels = "16 Data Points", xlabel = L"t", ylabel = "Incidence (Vector)",dpi = 400, guidefontfamily = "Computer Modern")
plot!(times, sol30Vector, legend=:topright, labels = "31 Data Points", xlabel = L"t", dpi = 400)
plot!(times, solTrueVector, legend=:topright, labels = "True Parameters", xlabel = L"t", dpi = 400)
scatter!(times, noisyDataVector, labels = "Noisy Incidence Vector Data")
display(plt2)
# savefig(plt2, "vectorIncidenceFittedTrue")

R_0(beta_v1, beta_h1, mu_h1, gamma1, mu_v1) = sqrt((beta_h1 * beta_v1 * 1000 * 5000)/((mu_h1 + gamma1) * mu_v1))

g(beta_v1, betah_1, gamma1) = R_0(beta_v1, betah_1, 0.00004, gamma1,0.09)

println("True: ", g(p0[2], p0[1], p0[3]))
println("16 data points: ", g(p15[2], p15[1], p15[3]))
println("31 data points: ", g(p30[2], p30[1], p30[3]))

trueR_0 = g(p0[2], p0[1], p0[3])
R_016 =  g(p15[2], p15[1], p15[3])
R_031 = g(p30[2], p30[1], p30[3])

println("16 Data Points Relative Error: ", (R_016 - trueR_0)/(trueR_0))
println("31 Data Points Relative Error: ", (R_031 - trueR_0)/(trueR_0))