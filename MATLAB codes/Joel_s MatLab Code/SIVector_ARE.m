%% Testing Relative Error
% y = [1 2 ;10 0.1; 0.1 0.2; 0.3 0.06];
% y0 = [1 2; 9 0.09; 0.1 0.2; 0.3 0.05];
% RelativeError(y0,y)
%% Testing ARE
params = [1 2 3 4 5];
param_results = [1 2 3 4 6;1 2 3 4 6];
rel_error = abs(params(1:5) - param_results)./params(1:5)
avgrelerror = 100*mean(rel_error)
%% Testing
tspan = linspace(1,100,1000);
ARE([0.0001 0.1 0.00025 0.01 0.1],[0.0001 0.1 0.01 0.01 0.1],50,tspan,[0 0.01 0.05 0.1 0.2],[100 1000])
%% Functions
function error = RelativeError(y0,y)
    diffsq = (y - y0).^2;
    weight = 1./(y0.^2);
    error = sum(sum((diffsq.*weight)));
end
function err = Model_Error(y,tspan,params,pops)
    params = max(params,0);  %Ensures nonnegative parameters
    %Uses fixed values for pih, piv
    y0 = [pops(1)*.999 pops(1)*.001 pops(2)*.999 pops(2)*.001];%Initial Condition
    pih = pops(1)*params(1);
    piv = pops(2)*params(2);
    [~,y2] = ode45(@SIVector_Model,tspan,y0,[],[params pih piv]); %Solves ode
    err = RelativeError(y,y2(:,[2,4])); %Finds error between y and params with model
end
function param_results = param_estimates(sigma,prevalence,init_guess,num_iters,tspan,pops)
    param_results = zeros(num_iters,5);%Number of parameters to esimate
    parfor i=1:num_iters
        y_error = prevalence + prevalence.*normrnd(0,sigma,size(prevalence));%Creates noisy data, with relative error
        newpars = fminsearch(@(param)(Model_Error(y_error,tspan,param,pops)),init_guess); %Finds params
        param_results(i,:) = newpars;
    end
end
%Note: We assume constant population size (at least for now)
%params is the 5 we are measureing(same with init guess)
function are = ARE(params,init_guess,num_iters,tspan,error_levels,pops)
    %Create error-free data
    pih = pops(1)*params(1);
    piv = pops(2)*params(2);
    y0 = [pops(1)*.999 pops(1)*.001 pops(2)*.999 pops(2)*.001]; %Do we change based on the eigenvector?
    [~,y] = ode45(@SIVector_Model,tspan,y0,[],[params pih piv]);
    I = y(:,[2 4]);%Prevalaence data
    num_levels = length(error_levels);
    are = zeros([num_levels 5]);%Number of parameters being varied
    for i = 1:num_levels
        param_results = param_estimates(error_levels(i),I,init_guess,num_iters,tspan,pops);
        rel_error = abs(params(1:5) - param_results)./params(1:5);
        avgrelerror = 100*mean(rel_error);
        are(i,:) = avgrelerror;
    end
end