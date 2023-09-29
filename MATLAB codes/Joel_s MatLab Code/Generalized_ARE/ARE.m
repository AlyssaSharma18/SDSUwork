function [AREs,all_param_est] = ARE(model,observations,parameters,y0,init_guess,tspan,num_iters,error_levels,obj_fun,fminfun,noisemaker,isStiff)
    % The measurment of the difference between two sets of observations
    if (nargin<9)||isempty(obj_fun)
        obj_fun = @OLS;
    end
    % The function which implements an optimization algorithm
    if (nargin<10)||isempty(fminfun)
        options = optimset('Display','none');
        fminfun = @(fun,init)(fminsearch(fun,init,options));
    end
    % The function which adds noise to our data
    if (nargin<11)||isempty(noisemaker)
        noisemaker = @relative_noise;
    end
    % Function which determines if we should run ode45 to get actual error
    if (nargin<12)||isempty(isStiff)
        isStiff = @(params) false;
    end
    rng(5);
    [~,y] = ode45(model,tspan,y0,[],parameters);%Correct Data
    data = observations(y);
    num_levels = length(error_levels);
    num_params = length(parameters);
    AREs = zeros([num_levels num_params]);
    all_param_est = zeros([num_levels num_iters num_params]);
    % For each error level
    for i = 1:num_levels
        %Estimate parameters num_iters times
        param_estimates = param_estimate(model,observations,error_levels(i),data,y0,init_guess,tspan,num_iters,num_params,obj_fun,fminfun,noisemaker,isStiff);
        % Calculate the average relative error
        rel_error = abs(parameters - param_estimates)./parameters;
        are = 100*mean(rel_error);
        %Store
        AREs(i,:) = are;
        all_param_est(i,:,:) = param_estimates;
    end
end
function param_results = param_estimate(model,observations,sigma,data,y0,init_guess,tspan,num_iters,num_params,obj_fun,fminfun,noisemaker,isStiff)
    param_results = zeros([num_iters num_params]);
    %Setting rng for reproducibility
    seed = round(1000*sigma);
    sc = parallel.pool.Constant(RandStream('threefry','Seed',seed));
    parfor i=1:num_iters %Uses parallel processing
        % Sets rng value
        stream = sc.Value;
        stream.Substream = i;
        prev = RandStream.setGlobalStream(stream);
        %Creates noisy data
        error_data = noisemaker(sigma,data);
        %Defines error between parameter estimate and noisy data
        error_fun = @(params)(Model_Error(model,observations,error_data,params,y0,tspan,obj_fun,data,isStiff));
        %Finds parameters which minimize this error
        newpars = fminfun(error_fun,init_guess);
        param_results(i,:) = newpars;
        RandStream.setGlobalStream(prev);
    end
end
function err = Model_Error(model,observations,error_data,params,y0,tspan,obj_fun,true_data,isStiff)
    %If parameters are negative or stiff, do not run ode45
    if any(params<0) || isStiff(params)
        err = 1e10; %Try inf?
    else
        [~,y_new] = ode45(model,tspan,y0,[],params);
        data_new = observations(y_new);
        err = obj_fun(error_data,data_new,true_data);
    end
end
function er = OLS(data1,data2,~)
    er = sum((data1 - data2).^2,'all');
end
function noisy = relative_noise(sigma,data)
noisy = data + data.*normrnd(0,sigma,size(data));
end