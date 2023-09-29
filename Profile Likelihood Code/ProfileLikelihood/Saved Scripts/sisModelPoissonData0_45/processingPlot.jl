# Activating ProfileLikelihood package 
using Pkg
Pkg.activate("ProfileLikelihood")

# Packages 
using Revise, ProfileLikelihood
using DiffEqBase, OrdinaryDiffEq, Plots, LikelihoodProfiler, DataFrames, Random, Distributions, CSV, LaTeXStrings, BenchmarkTools, Measures, BlackBoxOptim, BenchmarkTools, JLD2, Metaheuristics, StatsPlots

f = jldopen("sisIncidencePoisson0_45.jld2", "r")
theta3 = read(f, "theta3")
sol3 = read(f, "sol3")
threshold = read(f, "threshold")
PLconst = read(f, "PLconst")
thresholdP = read(f, "thresholdP")
close(f)

PLgamma = plot(theta3, [sol3, (x) -> (threshold + PLconst), (x) -> (thresholdP + PLconst)], xlabel = L"\gamma", ylabel = L"\chi^2_{\rm PL}", yformatter = :plain, legend=:topright, labels = [L"\chi^2_{\rm PL}" "Simultaneous Threshold" "Pointwise Threshold"], right_margin=5mm, dpi = 400)
scatter!([paramsFitted[3]], [loss + PLconst], color = "orange", labels = "Fitted Parameter")
display(PLgamma)
savefig(PLgamma, "PLgamma.png")