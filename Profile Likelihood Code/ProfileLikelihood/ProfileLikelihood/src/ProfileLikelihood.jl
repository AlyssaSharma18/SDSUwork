module ProfileLikelihood

using DiffEqBase, OrdinaryDiffEq, Plots, NLopt, Random, Distributions, BlackBoxOptim, LaTeXStrings, Measures, FiniteDifferences, Metaheuristics
export DiffEqBase, OrdinaryDiffEq, Plots, NLopt, Random, Distributions, BlackBoxOptim, LaTeXStrings, Measures, FiniteDifferences, Metaheuristics

# Generate Data
function generateData(index, seed, dist, prob, times, solver_opts; incidenceStatus = false)
    if incidenceStatus == false
        # Solve ODE
        sol = solve(
            prob,
            solver_opts[:alg],
            reltol=solver_opts[:reltol],
            abstol=solver_opts[:abstol],
            saveat=0.01,
            save_idxs=index
        )
        perfectData = [sol(t) for t in times]
        perfectData = abs.(perfectData)
    else
        perfectData = generateIncidenceData(index, prob, times, solver_opts)
        perfectData = abs.(perfectData)
    end
    # Set seed for reproducibility 
    Random.seed!(seed)
    # Create probability distributions to generate noisy data
    noisyDataGenerator = [dist(i) for i in perfectData]
    # Generate noisy data 
    noisyData = Vector{Float64}()
    for i in eachindex(noisyDataGenerator) 
        append!(noisyData, Distributions.rand(noisyDataGenerator[i], 1))
    end
    return perfectData, noisyData
end

function generateIncidenceData(index, prob, times, solver_opts)
    sol = solve(
        prob,
        solver_opts[:alg],
        reltol=solver_opts[:reltol],
        abstol=solver_opts[:abstol],
        saveat=0.01,
        save_idxs=index
    )
    cumulativeData = [sol(t) for t in times]
    perfectData = Vector{Float64}()
    for i in 1:length(cumulativeData)
        if i == 1
            append!(perfectData, 0)
        else 
            append!(perfectData, sol(times[i]) - sol(times[i-1]))
        end
    end
    return perfectData
end 

function generateIncidenceData(sol)
    incidenceData = Vector{Float64}()
    for i in 1:length(sol)
        if i == 1
            append!(incidenceData, 0)
        else 
            append!(incidenceData, sol[i] - sol[i-1])
        end
    end
    return incidenceData
end
##################################################
## Functions for different likelihood functions ##
##################################################

function relativeError(data, sol, noiseLevel)
    return (1 / noiseLevel^2) * sum(((abs.(sol) - data) ./ abs.(sol)) .^ 2) + 2.0 * sum(log.(abs.(sol)))
end

function poissonError(data, sol)
    return 2*(sum(abs.(sol)) - sum(data .* log.(abs.(sol))))
end

function constVarianceError(data, sol, sigma)
    return (1 / sigma^2) * sum((sol - data) .^ 2)
end

# Likelihood function
function likelihood(paramsCur, data, solObserved, prob, solver_opts, times, obj; paramIndex=0, paramEval=0)
    if paramIndex != 0
        paramsCurCopy = copy(paramsCur)
        insert!(paramsCurCopy, paramIndex, paramEval)
        probCur = remake(prob, p=paramsCurCopy)
    else
        probCur = remake(prob, p=paramsCur)
    end

    # solve odes
    sol = solve(
        probCur,
        solver_opts[:alg],
        reltol=solver_opts[:reltol],
        abstol=solver_opts[:abstol],
        saveat=times,
        save_idxs=solObserved
    )
    sol = sol[1, :]
    # Loss
    return obj(data, sol)
end

# Multiple dispatch for likelihood
function likelihood(paramsCur, data::Vector{Vector{T}}, solObserved, prob, solver_opts, times, objArr; incidenceObserved = [], paramIndex=0, paramEval=0) where {T<:Real}
    solObservedCopy = vcat(solObserved, incidenceObserved)
    if paramIndex != 0
        paramsCurCopy = copy(paramsCur)
        insert!(paramsCurCopy, paramIndex, paramEval)
        probCur = remake(prob, p=paramsCurCopy)
    else
        probCur = remake(prob, p=paramsCur)
    end

    # solve odes
    sol = solve(
        probCur,
        solver_opts[:alg],
        reltol=solver_opts[:reltol],
        abstol=solver_opts[:abstol],
        saveat=times,
        save_idxs=solObservedCopy
    )
    # loss
    loss = 0.0
    index = 0 
    for ind in eachindex(solObserved)
        loss += objArr[ind](data[ind], sol[ind,:])
        index += 1
    end
    for ind in eachindex(incidenceObserved)
        incidenceSol = generateIncidenceData(sol[ind + index,:])
        loss += objArr[ind + index](data[ind + index][2:end], incidenceSol[2:end]) # remove the first data point for cumulative data since it is always 0 and taking the log of 0 is -Inf  
    end
    return loss
end

# Estimate parameters
function estimateParams(p0, fitter_opts, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = [], paramIndex=0, paramEval=0, order = 4, maxRange = 1e-4, status = "local")
    if paramIndex == 0
        g = (paramsCur) -> likelihood(paramsCur, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved)
        searchRange = fitter_opts[:searchRange]
        resFirst = bboptimize(g, p0; Method=fitter_opts[:algFirst], SearchRange=searchRange, MaxTime=fitter_opts[:maxTimeFirst], TraceMode=fitter_opts[:traceMode], PopulationSize=1000)
        return best_fitness(resFirst), best_candidate(resFirst)
    else
        if status == "local"
        g = (paramsCur) -> likelihood(paramsCur, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved, paramIndex=paramIndex, paramEval=paramEval)
        function h(p, gradient)
            if length(gradient) > 0
                gradient[:] = grad(central_fdm(order, 1, max_range=maxRange), g, p)[1]
            end
            return g(p)
        end
        println("Gradient: ", grad(central_fdm(order, 1, max_range=maxRange), g, p0)[1])
        opt = Opt(fitter_opts[:algNLOpt], length(p0))
        opt.min_objective = (p, grad) -> h(p, grad)
        lowerBoundsCopy = copy(fitter_opts[:lowerBoundsNLOpt])
        upperBoundsCopy = copy(fitter_opts[:upperBoundsNLOpt])
        lowerBounds = deleteat!(lowerBoundsCopy, paramIndex)
        upperBounds = deleteat!(upperBoundsCopy, paramIndex)
        opt.lower_bounds = lowerBounds
        opt.upper_bounds = upperBounds
        opt.ftol_rel = 1e-7
        opt.ftol_abs = 1e-14
        (minf, minx, ret) = NLopt.optimize(opt, p0)
        println("Loss: ", minf)
        println("Fitted Parameters: ", minx)
        println("Return value: ", ret)
        # searchRangeCopy = copy(fitter_opts[:searchRange])
        # searchRange = deleteat!(searchRangeCopy, paramIndex)
        # resFirst = bboptimize(g, minx; Method=fitter_opts[:algFirst], SearchRange=searchRange, MaxTime=fitter_opts[:maxTimeFirst], TraceMode=fitter_opts[:traceMode], PopulationSize=500)
        # return best_fitness(resFirst), best_candidate(resFirst)
        return minf, minx
    elseif status == "global"
        g = (paramsCur) -> likelihood(paramsCur, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved, paramIndex=paramIndex, paramEval=paramEval)
        bounds = copy(fitter_opts[:boundsMH])
        bounds = bounds[:, 1:end .!= paramIndex]
        println(bounds)
        options = Options(time_limit = 240.0, f_calls_limit = 10000000000, iterations = 200000000)
        information = Information()
        algor = DE(N = 100, options = options, information = information)
        resOpt = Metaheuristics.optimize(g, bounds, algor)
        println(resOpt)
        return minimum(resOpt), minimizer(resOpt)
    elseif status == "globalBB"
        g = (paramsCur) -> likelihood(paramsCur, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved, paramIndex=paramIndex, paramEval=paramEval)
        searchRange = copy(fitter_opts[:searchRange])
        deleteat!(searchRange, paramIndex)
        resFirst = bboptimize(g, p0; Method=fitter_opts[:algFirst], SearchRange=searchRange, MaxTime=fitter_opts[:maxTimeFirst], TraceMode=fitter_opts[:traceMode], PopulationSize=1000)
        return best_fitness(resFirst), best_candidate(resFirst)
    end
    end
end

# Find threshold 
function findThresholdOptimization(confidence, numsParams, loss)
    threshold = loss + cquantile(Chisq(numsParams), 1 - confidence)
    println("The threshold to use for optimization is $threshold.")
    return threshold
end

##############################################
## Functions for finding profile likelihood ##
##############################################
function minPoint(paramIndex, loss, parametersFitted)
    return [parametersFitted[paramIndex]], [loss]
end

function goRightPL(stepSize, maxSteps, paramIndex, parametersFitted, data, solObserved, prob, solver_opts, times, loss, upperBound, fitter_opts, obj, incidenceObserved, order, maxRange, status)
    thetaRight = Vector{Float64}()
    solRight = Vector{Float64}()
    currVal = loss
    paramsFittedCopy = copy(parametersFitted)
    paramToLookAt = paramsFittedCopy[paramIndex]
    paramGuess = deleteat!(paramsFittedCopy, paramIndex)
    iter = 0
    if paramToLookAt + stepSize > fitter_opts[:searchRange][paramIndex][2]
        return thetaRight, solRight
    end
    while currVal < upperBound && iter < maxSteps && paramToLookAt < fitter_opts[:searchRange][paramIndex][2]
        paramToLookAt = paramToLookAt + stepSize
        append!(thetaRight, paramToLookAt)
        loss, paramGuess = estimateParams(paramGuess, fitter_opts, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved, paramIndex=paramIndex, paramEval=paramToLookAt, order = order, maxRange = maxRange, status = status)
        append!(solRight, loss)
        iter = iter + 1
        currVal = loss
    end
    return thetaRight, solRight
end

function goLeftPL(stepSize, maxSteps, paramIndex, parametersFitted, data, solObserved, prob, solver_opts, times, loss, upperBound, fitter_opts, obj, incidenceObserved, order, maxRange, status)
    thetaLeft = Vector{Float64}()
    solLeft = Vector{Float64}()
    currVal = loss
    paramsFittedCopy = copy(parametersFitted)
    paramToLookAt = paramsFittedCopy[paramIndex]
    paramGuess = deleteat!(paramsFittedCopy, paramIndex)
    iter = 0
    if paramToLookAt - stepSize < fitter_opts[:searchRange][paramIndex][1]
        return reverse(thetaLeft), reverse(solLeft)
    end
    while currVal < upperBound && iter < maxSteps && paramToLookAt > fitter_opts[:searchRange][paramIndex][1]
        paramToLookAt = paramToLookAt - stepSize
        append!(thetaLeft, paramToLookAt)
        loss, paramGuess = estimateParams(paramGuess, fitter_opts, data, solObserved, prob, solver_opts, times, obj; incidenceObserved = incidenceObserved, paramIndex=paramIndex, paramEval=paramToLookAt, order = order, maxRange = maxRange, status = status)
        append!(solLeft, loss)
        iter = iter + 1
        currVal = loss
    end
    return reverse(thetaLeft), reverse(solLeft)
end

function PL(stepSize, maxSteps, paramIndex, parametersFitted, data, solObserved, upperBound, loss, prob, solver_opts, times, fitter_opts, obj; incidenceObserved = [], order = 4, maxRange = 1e-4, status = "local")
    thetaRight, solRight = goRightPL(stepSize, maxSteps, paramIndex, parametersFitted, data, solObserved, prob, solver_opts, times, loss, upperBound, fitter_opts, obj, incidenceObserved, order, maxRange, status)
    thetaLeft, solLeft = goLeftPL(stepSize, maxSteps, paramIndex, parametersFitted, data, solObserved, prob, solver_opts, times, loss, upperBound, fitter_opts, obj, incidenceObserved, order, maxRange, status)
    thetaPoint, solPoint = minPoint(paramIndex, loss, parametersFitted)
    theta = vcat(thetaLeft, thetaPoint, thetaRight)
    sol = vcat(solLeft, solPoint, solRight)
    return theta, sol
end

function likelihoodConst(obj; noiseLevel=0.01, times=Vector{Float64}(), data=Vector{Float64}(), sigma = 1)
    if obj == "relativeError"
        return length(times) * log(1 / noiseLevel^2) + length(times) * log(2 * pi)
    elseif obj == "poissonError"
        return 2*sum(log.(factorial.(big.(data))))
    elseif obj == "constVarianceError"
        return length(times) * log(2 * pi) + length(times) * log(sigma^2)
    else
        return 0.0
    end
end

end # module 