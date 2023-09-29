# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, MATLAB

f = jldopen("testingDataCollection.jld2", "r")
bbLossArr = read(f, "bbLossArr")
trueLossArr = read(f, "trueLossArr")
mhLossArr = read(f, "mhLossArr")
bbParamsArr = read(f, "bbParamsArr")
mhParamsArr = read(f, "mhParamsArr")
closeOkay = read(f,"close")
close(f)

overfitArr = trueLossArr - bbLossArr
println(length(overfitArr))
mat"""
    $u = mle($(overfitArr),'pdf',@(x,a)gampdf(x,a,2),'start',1.5)
"""
println(2 * u)