# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, MATLAB, StatsPlots

f = jldopen("collectionData0_15.jld2", "r")
bbLossArr = read(f, "bbLossArr")
trueLossArr = read(f, "trueLossArr")
bbParamsArr = read(f, "bbParamsArr")
close(f)

f1 = jldopen("testingDataCollection.jld2", "r")
bbLossArr1 = read(f1, "bbLossArr")
trueLossArr1 = read(f1, "trueLossArr")
bbParamsArr1 = read(f1, "bbParamsArr")
close(f1)

bbLossArr1 = 2 * bbLossArr1
trueLossArr1 = 2 * trueLossArr1

append!(bbLossArr, bbLossArr1)
append!(trueLossArr, trueLossArr1)
overfitArr = trueLossArr - bbLossArr
println("The length of the array is ", length(overfitArr))
println(bbParamsArr)
# overfitArr = overfitArr[1:75]
println(length(overfitArr))
mat"""
    [$u, $d] = mle($(overfitArr),'pdf',@(x,a)gampdf(x,a,2),'start',1.5)
"""
println(2 * u)
println(2 .* d)
