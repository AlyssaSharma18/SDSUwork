%% Changing Number of Observations
num_points = linspace(5,101,97);
are = zeros(length(num_points),3);
model = @fixedpimu;
prev = @prev_both;
params = [0.0001 0.001 0.08996];
init_guess = params;
y0 = [99 1 999 1];
num_iters = 10;
error_levels = 0.1;
%loop
for i = 1:length(num_points)
    tspan = linspace(0,1000,num_points(i));
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting changing number of observations
figure
hold on
plot(num_points,are(:,1),'LineWidth',1.5);
plot(num_points,are(:,2),'LineWidth',1.5);
plot(num_points,are(:,3),'LineWidth',1.5);
plot(num_points,10*ones(length(num_points)),'--');
xlabel('Number of Observations');
ylabel('Average Relative Error');
legend('\beta_h','\beta_v','\gamma');
title('ARE for Different Numbers of Observations at \sigma = 10%');
hold off
% saveas(gcf,'SIV_ARE_Num_Obs.png')
%% Changing end date
end_dates = linspace(100,1000,101);
are = zeros(length(end_dates),3);
model = @fixedpimu;
prev = @prev_both;
params = [0.0001 0.001 0.08996];
init_guess = params;
y0 = [99 1 999 1];
num_iters = 10;
error_levels = 0.1;
for i = 1:length(end_dates)
    tspan = linspace(0,end_dates(i),101);
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting change in end date
figure
% model = @fixedpimu;
tspan_precise = linspace(0,1000,100);
%[t_precise,y_precise] = ode45(model,tspan_precise,y0,[],params);
hold on
plot(end_dates,are(:,1),'LineWidth',1.5)
plot(end_dates,are(:,2),'LineWidth',1.5)
plot(end_dates,are(:,3),'LineWidth',1.5)
%plot(t_precise,y_precise(:,2),'LineWidth',1.5);
plot(tspan_precise,10*ones(length(tspan_precise)),'--')
xlabel('End Date');
ylabel('Average Relative Error');
legend('betah','betav','gamma');
title('ARE for Different End Dates at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_End_Date.png');
%% Changing Type of Observations
error_levels = linspace(0,.3,50);
are = zeros(3,length(error_levels),3);
prev_b = @prev_both;
prev_h = @prev_host;
prev_v = @prev_vec;
model = @fixedpimu;
prev = @prev_both;
params = [0.0001 0.001 0.08996];
init_guess = params;
y0 = [99 1 999 1];
tspan = linspace(0,1000,101);
num_iters = 10;
are(1,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(2,:,:) = ARE(model,prev_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(3,:,:) = ARE(model,prev_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
%% Plotting changing Type of observations
%Add incidence or cumulative later, perhaps
%beta_h
figure
hold on
plot(error_levels,are(1,:,1),'LineWidth',1.5);
plot(error_levels,are(2,:,1),'LineWidth',1.5);
plot(error_levels,are(3,:,1),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('Both Prevalence','Host Prevalence','Vector Prevalence');
title('ARE of \beta_h for Types of Observations');
hold off
% saveas(gcf,'SIV_ARE_Obs_type_betah.png');
%beta_v
figure
hold on
plot(error_levels,are(1,:,2),'LineWidth',1.5);
plot(error_levels,are(2,:,2),'LineWidth',1.5);
plot(error_levels,are(3,:,2),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('Both Prevalence','Host Prevalence','Vector Prevalence');
title('ARE of \beta_v for Types of Observations');
hold off
% saveas(gcf,'SIV_ARE_Obs_type_betav.png');
%gamma
figure
hold on
plot(error_levels,are(1,:,3),'LineWidth',1.5);
plot(error_levels,are(2,:,3),'LineWidth',1.5);
plot(error_levels,are(3,:,3),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('Both Prevalence','Host Prevalence','Vector Prevalence');
title('ARE of \gamma for Types of Observations');
hold off
% saveas(gcf,'SIV_ARE_Obs_type_gamma.png');
%% Changing Optimization Algorithm
error_levels = linspace(0,.3,50);
are = zeros(3,length(error_levels),3);
params = [0.0001 0.001 0.08996];
init_guess = params;
model = @fixedpimu;
prev = @prev_both;
y0 = [99 1 999 1];
tspan = linspace(0,1000,101);
num_iters = 10;
are(1,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]); %fminsearch
are(2,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],@fmcon);%fmincon
are(3,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],@fmunc);%fminsearchbnd breaks stuff
%% Plot changing optimization algorithm
%beta_h
figure
hold on
plot(error_levels,are(1,:,1),'LineWidth',1.5);
plot(error_levels,are(2,:,1),'LineWidth',1.5);
plot(error_levels,are(3,:,1),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('fminsearch','fmincon','fminunc');
title('ARE of \beta_h for Optimization Algorithms');
hold off
%saveas('SIV_ARE_Opt_Algo_betah.png')
%beta_v
figure
hold on
plot(error_levels,are(1,:,2),'LineWidth',1.5);
plot(error_levels,are(2,:,2),'LineWidth',1.5);
plot(error_levels,are(3,:,2),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('fminsearch','fmincon','fminunc');
title('ARE of \beta_v for Optimization Algorithms');
hold off
%saveas('SIV_ARE_Opt_Algo_betav.png')
%gamma
figure
hold on
plot(error_levels,are(1,:,3),'LineWidth',1.5);
plot(error_levels,are(2,:,3),'LineWidth',1.5);
plot(error_levels,are(3,:,3),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
legend('fminsearch','fmincon','fminunc');
title('ARE of \gamma for Optimization Algorithms');
hold off
%saveas('SIV_ARE_Opt_Algo_gamma.png')
%% Changing Objective Functions
error_levels = linspace(0,.3,5);
are = zeros(3,length(error_levels),3);
prev_b = @prev_both;
model = @fixedpimu;
prev = @prev_both;
params = [0.0001 0.001 0.08996];
init_guess = params;
y0 = [99 1 999 1];
tspan = linspace(0,1000,101);
num_iters = 10;
ols = @(data1,data2)(WLS(data1,data2,1,0));
poisson = @(data1,data2)(WLS(data1,data2,1,1/2));
gls = @(data1,data2)(WLS(data1,data2,1,1));
are(1,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,ols);
are(2,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,poisson);
are(3,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,gls);
%% Plot changing objective function
%beta_h
figure
hold on
plot(error_levels,are(1,:,1),'LineWidth',1.5);
plot(error_levels,are(2,:,1),'LineWidth',1.5);
plot(error_levels,are(3,:,1),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
xlabel('Error Level');
ylabel('Average Relative Error');
ylim([0 50]);
legend('OLS','Poisson','GLS');
title('ARE of \beta_h for Objective Functions');
hold off
%saveas('SIV_ARE_Obj_Fun_betah.png')
%beta_v
figure
hold on
plot(error_levels,are(1,:,2),'LineWidth',1.5);
plot(error_levels,are(2,:,2),'LineWidth',1.5);
plot(error_levels,are(3,:,2),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
ylim([0 50]);
xlabel('Error Level');
ylabel('Average Relative Error');
legend('OLS','Poisson','GLS');
title('ARE of \beta_v for Objective Functions');
hold off
%saveas('SIV_ARE_Obj_Fun_betav.png')
%gamma
figure
hold on
plot(error_levels,are(1,:,3),'LineWidth',1.5);
plot(error_levels,are(2,:,3),'LineWidth',1.5);
plot(error_levels,are(3,:,3),'LineWidth',1.5);
plot(error_levels,100*error_levels,'--');
ylim([0 50]);
xlabel('Error Level');
ylabel('Average Relative Error');
legend('OLS','Poisson','GLS');
title('ARE of \gamma for Objective Functions');
hold off
%saveas('SIV_ARE_Obj_Fun_gamma.png')
%% Functions
function data = prev_host(y)
    data = y(:,2);
end
function data = prev_vec(y)
    data = y(:,4);
end
function data = prev_both(y)
data = y(:,[2 4]);
end
function noisy_data = Weight_Error(sigma,data,weight,p)
    noisy_data = data + weight.*(data.^p).*normrnd(0,sigma,size(data));
end
function fval = WLS(data1,data2,weight,p)
    fval = sum(weight.^(-2).*(data1.^(-2*p)).*(data1-data2).^2,'all');
end
function newpars = fmcon(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);%For some reason, giving it 0 as a lower bound makes it significantly worse
end
function newpars = fmunc(error_fun,init_guess)
    options = optimoptions('fminunc','Display','none');
    newpars = fminunc(error_fun,init_guess,options);
end
function dx = fixedpimu(t,x,z)
    dx = SIVector_Model(t,x,[0.00004 0.09 z(1) z(2) z(3) 0.004 90]);
end