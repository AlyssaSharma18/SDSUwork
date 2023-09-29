%Using fminsearch and OLS, with relative error
%% Basic Plot
model = @TuncerSIR_model;
params = [0.25 0.5];
y0 = [.99 .01 0];
tspan = linspace(0,50,51);
tspan_precise = linspace(0,50,1000);
[t,y] = ode45(model,tspan,y0,[],params);
figure
hold on
plot(t,y(:,2),'o');
[t_precise,y_precise] = ode45(model,tspan_precise,y0,[],params);
plot(t_precise,y_precise(:,2));
xlabel('Day');
ylabel('Prevalence');
legend('data points');
title('Prevalence Over Time')
hold off
% saveas(gcf,'SIR_Prev_Over_Time.png')
%% Visualizing for fewer points
model = @SIR_model;
params = [0.25 0.5];
y0 = [.99 0.01 0];
tspan = linspace(0,50,3);
tspan_precise = linspace(0,50,1000);
[t,y] = ode45(model,tspan,y0,[],params);
figure
hold on
plot(t,y(:,2),'o');
[t_precise,y_precise] = ode45(model,tspan_precise,y0,[],params);
plot(t_precise,y_precise(:,2));
xlabel('Day');
ylabel('Prevalence');
legend('data points');
title('Prevalence Over Time')
hold off
% saveas(gcf,'SIR_Prev_Over_Time')
%% Changing Number of Observations
%Can always consider even fewer points, etc
num_points = linspace(3,51,49);
are = zeros(length(num_points),2);
for i = 1:length(num_points)
    num = num_points(i);
    model = @SIR_model;
    prev = @prevalence;
    params = [0.25 0.5];
    y0 = [.99 0.01 0];
    init_guess = [0.1 0.1];
    tspan = linspace(0,50,num);
    error_levels = 0.1;
    num_iters = 10;
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels);
end
%% Plot Changing number of points
%Seems like it is doing too good
figure
hold on
plot(num_points,are(:,1),'LineWidth',1.5)
plot(num_points,are(:,2),'LineWidth',1.5)
plot(num_points,10*ones(length(num_points)),'--')
xlabel('Number of Observations');
ylabel('Average Relative Error');
legend('alpha','beta');
title('ARE for Different Numbers of Observations at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_Num_Obs.png');
%% Changing End Date
end_dates = linspace(10,51,42);
are = zeros(length(end_dates),2);
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
init_guess = [0.1 0.1];
error_levels = 0.1;
num_iters = 10;
for i = 1:length(end_dates)
    tspan = linspace(0,end_dates(i),51);
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting Changing End Date
%Note: Number of Points Remains the Same
figure
model = @SIR_model;
tspan_precise = linspace(0,50,1000);
[t_precise,y_precise] = ode45(model,tspan_precise,y0,[],params);
hold on
plot(end_dates,are(:,1),'LineWidth',1.5)
plot(end_dates,are(:,2),'LineWidth',1.5)
yline(10,'--');
ylabel('Average Relative Error');
yyaxis right
plot(t_precise,100.*y_precise(:,2),'LineWidth',1.5);
xlabel('End Date');
legend('alpha','beta','Prevalence(100x)');
title('ARE for Different End Dates at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_End_Date.png');
%% Changing Start Date
start_dates = linspace(1,48,48);
are = zeros(length(start_dates),2);
model = @SIR_model;
prev = @(y)(y(2:end,2));
params = [0.25 0.5];
y0 = [.99 0.01 0];
init_guess = [0.1 0.1];
error_levels = 0.1;
num_iters = 50;
for i = 1:length(start_dates)
    tspan = [0 linspace(start_dates(i),50,51)];
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting Changing Start Date
%Note: Number of Points Remains the Same
figure
model = @SIR_model;
tspan_precise = linspace(0,50,1000);
[t_precise,y_precise] = ode45(model,tspan_precise,y0,[],params);
hold on
plot(start_dates,are(:,1),'LineWidth',1.5)
plot(start_dates,are(:,2),'LineWidth',1.5)
yline(10,'--');
ylabel('Average Relative Error');
yyaxis right
plot(t_precise,100.*y_precise(:,2),'LineWidth',1.5);
xlabel('Start Date');
legend('alpha','beta','Prevalence');
title('ARE for Different Start Dates at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_Start_Date.png');
%% Changing Initial Guess of Alpha(Small)
alphas = linspace(0.0001,0.01,100);
are = zeros(length(alphas),2);
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
error_levels = 0.1;
num_iters = 500;
tspan = linspace(0,50,51);
for i = 1:length(alphas)
    init_guess = [alphas(i) 0.1];
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plot
figure
hold on
plot(alphas,are(:,1),'LineWidth',1.5)
plot(alphas,are(:,2),'LineWidth',1.5)
plot(alphas,10*ones(length(alphas)),'--')
xlabel('\alpha');
ylabel('Average Relative Error');
legend('alpha','beta');
title('ARE for Different \alpha values at \sigma = 10%');
hold off
saveas(gcf,'SIR_ARE_small_alphas.png');
%% Changing Initial Guess of Alpha(Large)
alphas = linspace(3,5,50);
are = zeros(length(alphas),2);
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
error_levels = 0.1;
num_iters = 100;
tspan = linspace(0,50,51);
for i = 1:length(alphas)
    init_guess = [alphas(i) 0.1];
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting
figure
hold on
plot(alphas,min(are(:,1),100),'LineWidth',1.5)
plot(alphas,min(are(:,2),100),'LineWidth',1.5)
plot(alphas,10*ones(length(alphas)),'--')
xlabel('\alpha');
ylabel('Average Relative Error');
legend('alpha','beta');
title('ARE for Different \alpha values at \sigma = 10%');
hold off
saveas(gcf,'SIR_ARE_large_alphas.png');
%% Changing Initial Guess of Beta(Small)
betas = linspace(0.0001,0.01,100);
are = zeros(length(betas),2);
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
error_levels = 0.1;
num_iters = 10;
tspan = linspace(0,50,51);
for i = 1:length(betas)
    init_guess = [0.1 betas(i)];
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting
figure
hold on
plot(betas,are(:,1),'LineWidth',1.5)
plot(betas,are(:,2),'LineWidth',1.5)
plot(betas,10*ones(length(betas)),'--')
xlabel('\beta');
ylabel('Average Relative Error');
legend('alpha','beta');
title('ARE for Different \beta values at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_small_betas.png');
%% Changing Initial Guess of Beta(Large)
betas = linspace(7,8,1000);
are = zeros(length(betas),2);
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
error_levels = 0.1;
num_iters = 10;
tspan = linspace(0,50,51);
for i = 1:length(betas)
    init_guess = [0.1 betas(i)];
    are(i,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
end
%% Plotting Large Betas
figure
hold on
plot(betas,min(are(:,1),100),'LineWidth',1.5)
plot(betas,min(are(:,2),100),'LineWidth',1.5)
plot(betas,10*ones(length(betas)),'--')
xlabel('\beta');
ylabel('Average Relative Error');
legend('alpha','beta');
title('ARE for Different \beta values at \sigma = 10%');
hold off
% saveas(gcf,'SIR_ARE_large_betas.png');
%% Type of Observation
error_levels = linspace(0,.3,50);
are = zeros(4,length(error_levels),2);
init_guess = [0.1 0.1];
model = @SIR_model;
prev = @prevalence;
inc = @incidence;
cumul = @cumulative;
params = [0.25 0.5];
y0 = [.99 0.01 0];
num_iters = 10;
tspan = linspace(0,50,51);
are(1,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(2,:,:) = ARE(model,inc,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(3,:,:) = ARE(model,cumul,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(4,:,:) = ARE(model,cumul,params,y0,init_guess,tspan,num_iters,error_levels,[],[],@err_in_incid);
%% Plotting Type of Observation
error_levels_p = 100*error_levels;
save('SIR_ARE_Obs_type.mat','ares');
figure
hold on
plot(error_levels_p,are(1,:,1),'LineWidth',1.5)
plot(error_levels_p,are(2,:,1),'LineWidth',1.5)
plot(error_levels_p,are(3,:,1),'LineWidth',1.5)
plot(error_levels_p,are(4,:,1),'LineWidth',1.5)
plot(error_levels_p,error_levels_p,'--')
xlabel('Error Level (%)');
ylabel('Average Relative Error');
legend('Prevalence','Incidence','Cumulative','Cumulative(Incidence Error)');
title('ARE of \alpha for Types of Observations');
hold off
% saveas(gcf,'SIR_ARE_Obs_Type_a.png')
figure
hold on
plot(error_levels_p,are(1,:,2),'LineWidth',1.5)
plot(error_levels_p,are(2,:,2),'LineWidth',1.5)
plot(error_levels_p,are(3,:,2),'LineWidth',1.5)
plot(error_levels_p,are(4,:,2),'LineWidth',1.5)
plot(error_levels_p,error_levels_p,'--')
xlabel('Error Level (%)');
ylabel('Average Relative Error');
legend('Prevalence','Incidence','Cumulative','Cumulative(Incidence Error)');
title('ARE of \beta for Types of Observations');
hold off
% saveas(gcf,'SIR_ARE_Obs_Type_b.png')
%% Error Function
error_levels = linspace(0,.3,10);
are = zeros(2,length(error_levels),2);
init_guess = [0.1 0.1];
model = @SIR_model;
prev = @prevalence;
inc = @incidence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
num_iters = 100;
tspan = linspace(0,50,51);
are(1,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
are(2,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,@relative_error,[]);
%% Plotting
%Perhaps add the line
figure
hold on
plot(error_levels,are(1,:,1),'LineWidth',1.5)
plot(error_levels,are(1,:,2),'LineWidth',1.5)
plot(error_levels,are(2,:,1),'LineWidth',1.5)
plot(error_levels,are(2,:,2),'LineWidth',1.5)
plot(error_levels,100*error_levels,'--')
xlabel('Error Level');
ylabel('Average Relative Error(%)');
legend('OLS:\alpha','OLS:\beta','GLS:\alpha','GLS:\beta');
title('OLS vs GLS');
hold off
%saveas(gcf,'SIR_ARE_Obj_Fun.png')
%% FMinFun
error_levels = linspace(0,.3,50);
are = zeros(3,length(error_levels),2);
init_guess = [0.1 0.1];
model = @SIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = [.99 0.01 0];
num_iters = 20;
tspan = linspace(0,50,51);
are(1,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]); %fminsearch
are(2,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],@fmcon);%fmincon
are(3,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],@fminunc);%fminsearchbnd breaks stuff
% are(4,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],@fminsearchbnd);

%% Plotting
%Could change to % for error levels in future
figure
hold on
plot(error_levels,are(1,:,1),'b','LineWidth',1.5)
plot(error_levels,are(2,:,1),'r','LineWidth',1.5)
plot(error_levels,are(3,:,1),'m','LineWidth',1.5)
plot(error_levels,are(1,:,2),'b--','LineWidth',1.5)
plot(error_levels,are(2,:,2),'r--','LineWidth',1.5)
plot(error_levels,are(3,:,2),'m--','LineWidth',1.5)
xlabel('Error Level(%)');
ylabel('Average Relative Error(%)');
legend('fminsearch:\alpha','fmincon:\alpha','fminunc:\alpha','fminsearch:\beta','fmincon:\beta','fminunc:\beta');
title('ARE for Different Optimizers');
hold off
%saveas(gcf,'SIR_ARE_Error_Fun.png')
%% Plotting UNUSED
% %Could change to % for error levels in future
% figure
% hold on
% plot(error_levels,are(1,:,1),'LineWidth',1.5)
% plot(error_levels,are(2,:,1),'LineWidth',1.5)
% plot(error_levels,are(3,:,1),'LineWidth',1.5)
% plot(error_levels,are(4,:,1),'LineWidth',1.5)
% xlabel('Error Level(%)');
% ylabel('Average Relative Error(%)');
% legend('fminsearch','fmincon','fminunc','fminsearchbnd');
% title('ARE of \alpha for Different Optimizers');
% hold off
% figure
% hold on
% plot(error_levels,are(1,:,2),'LineWidth',1.5)
% plot(error_levels,are(2,:,2),'LineWidth',1.5)
% plot(error_levels,are(3,:,2),'LineWidth',1.5)
% plot(error_levels,are(4,:,2),'LineWidth',1.5)
% xlabel('Error Level(%)');
% ylabel('Average Relative Error(%)');
% legend('fminsearch','fmincon','fminunc','fminsearchbnd');
% title('ARE of \beta for Different Optimizers');
% hold off
% %saveas(gcf,'SIR_ARE_Error_Fun.png')
%% Testing
y = [90 10 0; 89 11 0; 87 12 1; 80 18 2;];
c = cumulative(y)
inc = incidence(y)
cumul_to_incid(c)
incid_to_cumul(inc)
err_in_incid(0.1,cumulative(y))
%% Functions
function data = prevalence(y)
    data = y(:,2);
end
% function data = cont_incidence_UNUSED(y)
%  data = y(:,1).*y(:,2).*0.5;
% end
function data = cumulative(y)
    start = y(1,2) + y(1,3);
    data = y(:,2)+y(:,3)-start;
end
function data = incidence(y)
num_days = length(y(:,1));
data = zeros(num_days,1);
data(1) = 0;
for i = 2:length(y(:,1))
    data(i) = y(i,2) + y(i,3)-y(i-1,2)-y(i-1,3);
end
end
function data = cumul_to_incid(y)
    num_days = length(y);
    data = zeros(num_days,1);
    data(1) = 0;
    for i = 2:num_days
        data(i) = y(i)-y(i-1);
    end
end
function data = incid_to_cumul(y)
    num_days = length(y);
    data = zeros(num_days,1);
    data(1) = 0;
    for i = 2:num_days
        data(i) = y(i) + data(i-1);
    end
end
function noisy_data = err_in_incid(sigma,y)
    noisy_data = zeros(size(y));
    for i = 1:size(y,2)
        incid = cumul_to_incid(y(:,i));
        noisy_incid = incid + incid.*normrnd(0,sigma,size(incid));
        noisy_data(:,i) = incid_to_cumul(noisy_incid);
    end
end
function er = relative_error(data1,data2)
    diffsq = (data1 - data2).^2;
    weight = 1./(data1.^2);
    er = sum(sum((diffsq.*weight)));%Use dimensions for preciseness-this is treated as a vector, and a norm
end
function newpars = fmcon(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);%For some reason, giving it 0 as a lower bound makes it significantly worse
end
