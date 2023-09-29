# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, MATLAB, StatsPlots

f = jldopen("collectionData0_45_relativeError.jld2", "r")
bbLossArr = read(f, "bbLossArr")
trueLossArr = read(f, "trueLossArr")
bbParamsArr = read(f, "bbParamsArr")
close(f)

overfitArr = trueLossArr - bbLossArr
println(length(overfitArr))
mat"""
    [$u, $d] = mle($(overfitArr),'pdf',@(x,a)gampdf(x,a,2),'start',1.5)
"""
println(2 * u)
println(2 .* d)

plt = histogram(overfitArr, normalize = :probability, bins = 15, xlabel = L"\chi^2(\theta^{\ast}) - \chi^2(\theta)", ylabel = "proportion", fontfamily = "Computer Modern", labels = "Overfitted Data", dpi = 400)
plot!(Chisq(3.008147048950195), labels = L"$\chi^2_{3.001}$-distribution")
savefig(plt, "overfitPlot2.png")
[2.60198328231417, 3.4143108155862203]
# CI: 