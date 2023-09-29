%% Recovering Tuncer Results
%Parameters
alpha = 0.5;
beta = 1;
end_date = 100;
tspan = linspace(1,end_date,1000);
y0 = [.999 .001 0];
num_points = 100;
params = [alpha,beta];
[t,y] = ode45(@(t,x) SIR_model(t,x,params),tspan,y0);
figure
plot(t,y(:,2));
num_iter = 10;
ares = ARE_p([alpha,beta],[0.1 0.1],num_iter,num_points,end_date,[0 0.01 0.05 0.1 0.2 0.3])
%% ARE with fixed beta = 0.5
beta = 0.5;
alphas = linspace(0.01,0.55,50);
num_alphas = length(alphas);
ares = zeros([2,num_alphas]);
max_inf_times = zeros(num_alphas,1);
max_inf = zeros(num_alphas,1);
end_date = 100;
num_points = 100;
tspan = linspace(1,end_date,1000);
y0 = [.999 .001 0];
for i = 1:num_alphas
    ares(:,i) = ARE_p([alphas(i),beta],100,num_points,0.1,end_date);
    [t,y] = ode45(@SIR_model,tspan,y0,[],[alphas(i),beta]);
    [max_inf(i),max_ind] = max(y(:,2));
    max_inf_times(i) = tspan(max_ind);
end
%% Data For Making a Contour Plot
[alpha,beta] = meshgrid(0.05:0.01:0.1,0.05:0.02:0.5);
s = size(alpha);
y = zeros([s 2]);
cols = s(2);
rows = s(1);
for c = 1:cols
    for r = 1:rows
        y(r,c,:) = ARE_p([alpha(1,c),beta(r,1)],100,50,[0.1],365);
    end
end
%% Plotting
figure
plot(alphas(1:50),ares(1,1:50))
title('are vs alpha')
figure
plot(alphas,max_inf_times)
title('peak infection time vs alpha')
figure
plot(alphas,max_inf)
title('peak infected vs alpha')
%% Plotting
figure
contourf(alpha,beta,y(:,:,1))
title('ARE of \alpha at \sigma = 10%');
xlabel('\alpha');
ylabel('\beta');
colorbar;
%% Plotting beta
figure
contourf(alpha,beta,y(:,:,1))
title('ARE of \beta at \sigma = 10%');
xlabel('\alpha');
ylabel('\beta');
colorbar;
%% Functions
function error = OLS(y,y2)
    diffsq = (y - y2).^2;
    error = sum(diffsq);
end
function err = Model_Error(y,num_points,end_date,params)
    if (params(1) < 0) || (params(2) < 0)
        err = 1000;
    else
    params = max(params,0);  %Ensures nonnegative parameters
    tspan = linspace(1,end_date,num_points);
    y0 = [.999 .001 0]; %Initial condition
    [t,y2] = ode45(@SIR_model,tspan,y0,[],params); %Solves ode
    err = OLS(y,y2(:,2)); %Finds error between y and params with model
    end
end
function [param_results,fvals] = param_estimate(sigma,prevalence,org_params,init_guess,num_iters,num_points,end_date)
    param_results = zeros(num_iters,2);
    fvals = zeros(num_iters,1);
    parfor i=1:num_iters
        y_error = prevalence + prevalence.*normrnd(0,sigma,size(prevalence));%Creates noisy data
        %init_guess = abs(normrnd(org_params,0.1*org_params))
        [newpars,fval] = fminsearch(@(param)(Model_Error(y_error,num_points,end_date,param)),init_guess);%Finds params
%         disp(newpars)
        param_results(i,:) = newpars;
        fvals(i) = fval;
    end
end
function are = ARE(param_results,org_params)
    rel_error = abs(org_params - param_results)./org_params;
    are = 100*mean(rel_error);
end
function are = ARE_p(params,init_guess,num_iters,num_points,end_date,error_levels)
    %Makes error-free data
    tspan = linspace(1,end_date,num_points);
    y0 = [.999 .001 0];
    [t,y] = ode45(@SIR_model,tspan,y0,[],params);
    %Prevalence data
    prev = y(:,2);
    num_levels = length(error_levels);
    AREs = zeros([num_levels 2]);
    %Caluclates AREs
    for i = 1:num_levels
        [param_results, fvals] = param_estimate(error_levels(i),prev,params,init_guess,num_iters,num_points,end_date);
        are = ARE(param_results,params);
        AREs(i,:) = are;
    end
    are = AREs;
end