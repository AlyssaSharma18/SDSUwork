# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, Roots, Interpolations

f = jldopen("sisIncidencePoisson0_15.jld2", "r")
loss = read(f, "loss")
paramsFitted = read(f, "paramsFitted")
threshold = read(f, "threshold")
thresholdP = read(f, "thresholdP")
PLconst = read(f, "PLconst")
theta1 = read(f,"theta1")
sol1 = read(f, "sol1")
theta2 = read(f,"theta2")
sol2 = read(f, "sol2")
theta3 = read(f,"theta3")
sol3 = read(f, "sol3")
close(f)

# Interpolations 
interp_linear1 = LinearInterpolation(theta1, sol1)
interp_linear2 = LinearInterpolation(theta2, sol2)
interp_linear3 = LinearInterpolation(theta3, sol3)


# Threshold and Simultaneous Threshold
simThreshold = threshold + PLconst 
poiThreshold = thresholdP + PLconst

# Root finding 

# THETA_1 AND SOL_1
# Pointwise 
res1 = find_zeros((x) -> interp_linear1(x) - poiThreshold, (theta1[1], theta1[end]))
println("The pointwise confidence interval for sol1 is ", res1)

# Simultaneous
res2 = find_zeros((x) -> interp_linear1(x) - simThreshold, (theta1[1], theta1[end]))
println("The simultaneous confidence interval for sol1 is ", res2)

# THETA_2 AND SOL_2
# Pointwise 
res3 = find_zeros((x) -> interp_linear2(x) - poiThreshold, (theta2[1], theta2[end]))
println("The pointwise confidence interval for sol2 is ", res3)

# Simultaneous 
res4 = find_zeros((x) -> interp_linear2(x) - simThreshold, (theta2[1], theta2[end]))
println("The simultaneous confidence interval for sol2 is ", res4)

# THETA_3 AND SOL_3
# Pointwise 
res5 = find_zeros((x) -> interp_linear3(x) - poiThreshold, (theta3[1], theta3[end]))
println("The pointwise confidence interval for sol3 is ", res5)

# Simultaneous 
res6 = find_zeros((x) -> interp_linear3(x) - simThreshold, (theta3[1], theta3[end]))
println("The simultaneous confidence interval for sol3 is ", res6)