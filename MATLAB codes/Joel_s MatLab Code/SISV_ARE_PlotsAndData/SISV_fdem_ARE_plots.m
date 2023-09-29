%Do we label ARE as (%) in plots?
%% Setup
model = @SISV_fdem;
prev_b = @prev_both;
prev_h = @prev_host;
prev_v = @prev_vector;
muh = 0.00004;
muv = 0.1;
pih = muh;
piv = muv;
beta_h = 0.1;
beta_v = 0.2;
gamma = 0.09996;
params = [beta_h beta_v gamma];
% init_guess = params; %Always have correct initial guess
init_guess = [0.05 0.1 0.05];
y0 = [.99 .01 .99 0.01];
num_iters = 50;
error_levels = 0.1;
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
%% Optional Plotting stuff
r0 = sqrt(beta_h*beta_v/((0.00004+gamma)*0.1));
h_equ = 1-(muh + gamma + (beta_h/r0^2))/(muh + gamma + beta_h);
v_equ = 1-(muv + (beta_v/r0^2))/(muv+beta_v);
fprintf('R0: %f\nHost Equilibrium: %f\n Vector Equilibrium: %f\n',[r0 h_equ v_equ]);
[~,y] = ode45(model,tspan,y0,[],params);
figure
hold on
plot(tspan,y(:,[2 4]));
legend('Infected Hosts','Infected Vectors');
hold off
%% Changing number of Observations
num_obs = linspace(3,100,98);
ares = zeros(length(num_obs),3);
for i = 1:length(num_obs)
    tspan = linspace(0,e_date,num_obs(i));
    ares(i,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
end
%% Plotting Changing Number of Observations
figure
hold on
plot(num_obs,ares(:,1),'LineWidth',1.5);
plot(num_obs,ares(:,2),'LineWidth',1.5);
plot(num_obs,ares(:,3),'LineWidth',1.5);
yline(10,'--');
legend('\beta_h','\beta_v','\gamma');
xlabel('Number of Observations');
ylabel('Average Relative Error');
title('ARE for Number of Observations at \sigma=10%')
hold off
% resultsdata = [num_obs' ares];
% writematrix(resultsdata,'SISV_fdem_ARE_numobs.xlsx')
% saveas(gcf,'SISV_fdem_ARE_numobs.png');
%% Changing Final Date (1000 Observations)
num_obs = 1000;
final_dates = linspace(10,e_date,e_date-9);
ares = zeros(length(final_dates),3);
init_guess = [0.05 0.05 0.05];
pes = zeros(length(final_dates),num_iters,3);
for i = 1:length(final_dates)
    tspan = linspace(0,final_dates(i),num_obs);
    [ares(i,:),pes(i,:,:)] = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
end
%% Plotting Changing Final Date
%check double axis
[t,y] = ode45(@SISV_fdem,linspace(0,e_date,e_date+1),y0,[],params);
figure
hold on
plot(final_dates,ares(:,1),'LineWidth',1.5);
plot(final_dates,ares(:,2),'LineWidth',1.5);
plot(final_dates,ares(:,3),'LineWidth',1.5);
yline(10,'--');
yyaxis right
plot(t,100*y(:,2),'LineWidth',1.25);
plot(t,100*y(:,4),'LineWidth',1.25);
legend('\beta_h','\beta_v','\gamma','Host Prevalence(%)','Vector Prevalence(%)');
xlabel('Final Date of Observations');
ylabel('Average Relative Error');
title('ARE for Final Date of Observations at \sigma=10%');
hold off
% resultsdata = [final_dates' ares];
% writematrix(resultsdata,'SISV_fdem_ARE_findate.xlsx');
% saveas(gcf,'SISV_fdem_ARE_findate.xlsx');
%% Param Estimate Plotting
fig1 = Plot_Param_Estimates([pes(1,:,1);pes(1,:,2)],model,prev_b,[0.1 0.2],y0,tspan,@(p)([p(1) p(2) params(3)]),[]);
xlabel('\beta_h');ylabel('\beta_v');title('Parameter Estimates for Final Date of Day 10');
saveas(fig1,'SISV_PE_Final_Day_10_bhbv.png')
fig2 = Plot_Param_Estimates([pes(1,:,1);pes(1,:,3)],model,prev_b,[0.1 0.1],y0,tspan,@(p)([p(1) params(2) p(3)]),[]);
xlabel('\beta_h');ylabel('\gamma');title('Parameter Estimates for Final Date of Day 10');
saveas(fig2,'SISV_PE_Final_Day_10_bhg.png');
fig3 = Plot_Param_Estimates([pes(1,:,2);pes(1,:,3)],model,prev_b,[0.1 0.1],y0,tspan,@(p)([params(1) p(2) p(3)]),[]);
xlabel('\beta_v');ylabel('\gamma');title('Parameter Estimates for Final Date of Day 10');
saveas(fig3,'SISV_PE_Final_Day_10_bvg.png');
%% Changing Start Date
num_obs = 1000;
prev_b_span = @(y)(y(1:end,[2 4]));
start_dates = linspace(1,e_date-10,e_date-10);
ares = zeros(length(start_dates),3);
pes = zeros(length(start_dates),num_iters,3);
init_guess = [0.05 0.1 0.05];
for i = 1:length(start_dates)
    tspan = [0 linspace(start_dates(i),e_date,num_obs)];
    [ares(i,:), pes(i,:,:)] = ARE(model,prev_b_span,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
end
%% Plotting Changing Start Date
%check double axis
[t,y] = ode45(@SISV_fdem,linspace(0,e_date,e_date+1),y0,[],params);
figure
hold on
bh = plot(start_dates,ares(:,1),'LineWidth',1.5);
bv = plot(start_dates,ares(:,2),'LineWidth',1.5);
gm = plot(start_dates,ares(:,3),'LineWidth',1.5);
yline(10,'--');
yyaxis right
ih = plot(t,100*y(:,2),'LineWidth',1.25);
iv = plot(t,100*y(:,4),'LineWidth',1.25);
legend([bh,bv,gm,ih,iv],'\beta_h','\beta_v','\gamma','Host Prevalence(%)','Vector Prevalence(%)');
xlabel('Start Date of Observations');
ylabel('Average Relative Error');
title('ARE for Start Date of Observations at \sigma=10%');
hold off
% resultsdata = [start_dates' ares];
% writematrix(resultsdata,'SISV_fdem_ARE_startdate.xlsx');
% saveas(gcf,'SISV_fdem_ARE_startdate.xlsx');
%% Param Estimate Plotting
fig1 = Plot_Param_Estimates([pes(1,:,1);pes(1,:,2)],model,prev_b,[0.1 0.2],y0,tspan,@(p)([p(1) p(2) params(3)]),[]);
xlabel('\beta_h');ylabel('\beta_v');title('Parameter Estimates for Start Date of Day 10');
saveas(fig1,'SISV_PE_Start_Day_10_bhbv.png')
fig2 = Plot_Param_Estimates([pes(1,:,1);pes(1,:,3)],model,prev_b,[0.1 0.1],y0,tspan,@(p)([p(1) params(2) p(3)]),[]);
xlabel('\beta_h');ylabel('\gamma');title('Parameter Estimates for Start Date of Day 10');
saveas(fig2,'SISV_PE_Start_Day_10_bhg.png');
fig3 = Plot_Param_Estimates([pes(1,:,2);pes(1,:,3)],model,prev_b,[0.2 0.1],y0,tspan,@(p)([params(1) p(2) p(3)]),[]);
xlabel('\beta_v');ylabel('\gamma');title('Parameter Estimates for Start Date of Day 10');
saveas(fig3,'SISV_PE_Start_Day_10_bvg.png');
%% Changing Species of Observation
error_levels = linspace(0,.3,31);
ares = zeros(3,length(error_levels),3);
ares(1,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(2,:,:) = ARE(model,prev_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(3,:,:) = ARE(model,prev_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
%% Plotting Changing Species of Observation
error_levels_p = 100*error_levels;
%beta_h
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Both','Hosts Only','Vectors Only');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Observations of Different Species')
% saveddata = [error_levels_p; ares(:,:,1)]';
% writematrix(saveddata,'SISV_fdem_ARE_species_betah.xlsx')
% saveas(gcf,'SISV_fdem_ARE_species_betah.png')
%beta_v
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Both','Hosts Only','Vectors Only');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Observations of Different Species')%Need a better term, probably
% saveddata = [error_levels_p; ares(:,:,2)]';
% writematrix(saveddata,'SISV_fdem_ARE_species_betav.xlsx')
% saveas(gcf,'SISV_fdem_ARE_species_betav.png')
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Both','Hosts Only','Vectors Only');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \gamma for Observations of Different Species')%Need a better term, probably
% saveddata = [error_levels_p; ares(:,:,3)]';
% writematrix(saveddata,'SISV_fdem_ARE_species_gamma.xlsx')
% saveas(gcf,'SISV_fdem_ARE_species_gamma.png')
%% Changing Types of Observation for Hosts
model = @SISV_C_fdem;
y0 = [.99 .01 .99 .01 0 0];
incid_h = @incid_host;
cumul_h = @cumul_host;
error_levels = linspace(0,.3,31);
ares = zeros(3,length(error_levels),3);
ares(1,:,:) = ARE(model,prev_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(2,:,:) = ARE(model,incid_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(3,:,:) = ARE(model,cumul_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(4,:,:) = ARE(model,cumul_h,params,y0,init_guess,tspan,num_iters,error_levels,[],[],@err_in_incid,@under_hundred);
%% Plotting Changing Types of Observation for Hosts
save('SISV_ARE_hobs.mat','ares');
error_levels_p = 100*error_levels;
%beta_h
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Types of Host Observations')
% saveas(gcf,'SISV_fdem_ARE_hobs_types_betah.png')
%beta_v
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Types of Host Observations')
% saveas(gcf,'SISV_fdem_ARE_hobs_types_betav.png')
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \gamma for Types of Host Observations')
% saveas(gcf,'SISV_fdem_ARE_hobs_types_gamma.png')
%% Changing Types of Observation for Vectors
model = @SISV_C_fdem;
y0 = [.99 .01 .99 .01 0 0];
incid_v = @incid_vector;
cumul_v = @cumul_vector;
error_levels = linspace(0,.3,31);
ares = zeros(4,length(error_levels),3);
ares(1,:,:) = ARE(model,prev_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(2,:,:) = ARE(model,incid_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(3,:,:) = ARE(model,cumul_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(4,:,:) = ARE(model,cumul_v,params,y0,init_guess,tspan,num_iters,error_levels,[],[],@err_in_incid,@under_hundred);
%% Plotting Changing Types of Observation for Vectors
save('SISV_ARE_vobs.mat','ares');
error_levels_p = 100*error_levels;
%beta_h
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Types of Vector Observations')
% saveas(gcf,'SISV_fdem_ARE_vobs_types_betah.png')
%beta_v
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Types of Vector Observations')
% saveas(gcf,'SISV_fdem_ARE_vobs_types_betav.png')
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \gamma for Types of Vector Observations','Cumulative (Incidence Error)')
% saveas(gcf,'SISV_fdem_ARE_vobs_types_gamma.png')
%% Changing Types of Observation for Both
model = @SISV_C_fdem;
y0 = [.99 .01 .99 .01 0 0];
incid_b = @incid_both;
cumul_b = @cumul_both;
c_ierr_nmaker = @err_in_incid;
error_levels = linspace(0,.3,31);
ares = zeros(4,length(error_levels),3);
ares(1,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(2,:,:) = ARE(model,incid_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(3,:,:) = ARE(model,cumul_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
ares(4,:,:) = ARE(model,cumul_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],@c_ierr_nmaker,@under_hundred);
%% Plotting Changing Types of Observation for Both
error_levels_p = 100*error_levels;
%beta_h
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Types of Observations')
% saveddata = [error_levels_p; ares(:,:,1)]';
% writematrix(saveddata,'SISV_fdem_ARE_bobs_types_betah.xlsx')
% saveas(gcf,'SISV_fdem_ARE_bobs_types_betah.png')
%beta_v
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Types of Observations')
% saveddata = [error_levels_p; ares(:,:,2)]';
% writematrix(saveddata,'SISV_fdem_ARE_bobs_types_betav.xlsx')
% saveas(gcf,'SISV_fdem_ARE_bobs_types_betav.png')
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(4,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Prevalence','Incidence','Cumulative','Cumulative (Incidence Error)');
xlabel('\sigma(%)');
ylabel('Average Relative Error');
title('ARE of \gamma for Types of Observations');
% saveddata = [error_levels_p; ares(:,:,3)]';
% writematrix(saveddata,'SISV_fdem_ARE_bobs_types_gamma.xlsx')
% saveas(gcf,'SISV_fdem_ARE_bobs_types_gamma.png')
%% Changing Optimization Algorithm
% Many a possibility-start with the great ones
error_levels = linspace(0,.3,31);
num_iters = 10;
opt_algos = {@fminsearch,@fmcon,@fmunc};
ares = zeros([length(opt_algos), length(error_levels), 3]);
for i = 1:length(opt_algos)
    fminfun = opt_algos(i);
    ares(i,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],fminfun,[],@under_hundred);
end
%% Plotting Changes in Optimization Algorithm
error_levels_p = 100*error_levels;
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Nelder-Mead','interior-point','quasi-newton');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Optimization Algorithms')
hold off
% saveddata = [error_levels_p; ares(:,:,1)]';
% writematrix(saveddata,'SISV_fdem_ARE_Optim_Algo_betah.xlsx')
% saveas(gcf,'SISV_fdem_ARE_Optim_Algo_betah.png');
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Nelder-Mead','interior-point','quasi-newton');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Optimization Algorithms')
hold off
% saveddata = [error_levels_p; ares(:,:,2)]';
% writematrix(saveddata,'SISV_fdem_ARE_Optim_Algo_betav.xlsx')
% saveas(gcf,'SISV_fdem_ARE_Optim_Algo_betav.png');
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('Nelder-Mead','interior-point','quasi-newton');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \gamma for Optimization Algorithms')
hold off
% saveddata = [error_levels_p; ares(:,:,3)]';
% writematrix(saveddata,'SISV_fdem_ARE_Optim_Algo_gamma.xlsx')
% saveas(gcf,'SISV_fdem_ARE_Optim_Algo_gamma.png');
%% Plotting
error_levels_p = 100*error_levels;
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('100','50','10');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_h for Bound on Interior Point')
hold off
% saveddata = [error_levels_p; ares(:,:,1)]';
% writematrix(saveddata,'SISV_fdem_ARE_bounds_betah.xlsx')
% saveas(gcf,'SISV_fdem_ARE_bounds_betah.png');
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('100','50','10');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \beta_v for Bound on Interior Point')
hold off
% saveddata = [error_levels_p; ares(:,:,2)]';
% writematrix(saveddata,'SISV_fdem_ARE_bounds_betav.xlsx')
% saveas(gcf,'SISV_fdem_ARE_bounds_betav.png');
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
plot(error_levels_p,error_levels_p,'--');
legend('100','50','10');
xlabel('\sigma(%)')
ylabel('Average Relative Error')
title('ARE of \gamma for Bound on Interior Point')
hold off
% saveddata = [error_levels_p; ares(:,:,3)]';
% writematrix(saveddata,'SISV_fdem_ARE_Optim_Algo_gamma.xlsx')
% saveas(gcf,'SISV_fdem_ARE_Optim_Algo_gamma.png');
%% Observations Functions
function data = prev_both(y)
data = y(:,[2 4]);
end
function data = prev_host(y)
    data = y(:,2);
end
function data = prev_vector(y)
    data = y(:,4);
end
function data = cumul_both(y)
    data = y(:,[5 6]);
end
function data = cumul_host(y)
    data = y(:,5);
end
function data = cumul_vector(y)
    data = y(:,6);
end
function data = incid_host(y)
    num_days = length(y(:,1));
    data = zeros(num_days,1);
    data(1) = 0;
    for i = 2:length(y(:,1))
        data(i) = y(i,5)-y(i-1,5);
    end
end
function data = incid_vector(y)
    num_days = length(y(:,1));
    data = zeros(num_days,1);
    data(1) = 0;
    for i = 2:length(y(:,1))
        data(i) = y(i,6)-y(i-1,6);
    end
end
function data = incid_both(y)
    data = [(incid_host(y))' (incid_vector(y))'];
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
%% Optimization Functions
function newpars = fmcon(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);
end
function newpars = fmchundred(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[0 0 0],[100 100 100],[],options);
end
function newpars = fmcfifty(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[0 0 0],[50 50 50],[],options);
end
function newpars = fmcten(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[0 0 0],[10 10 10],[],options);
end
function newpars = fmunc(error_fun,init_guess)
    options = optimoptions('fminunc','Display','none');
    newpars = fminunc(error_fun,init_guess,options);
end
%% Setup Functions
function is = under_hundred(params)
    if (params(1) > 100 || params(2) > 100)
        is = true;
    else
        is = false;
    end
end
function date = end_date(model,y0,params)
    fulltspan = linspace(0,2000,10000);
    [t,y] = ode45(model,fulltspan,y0,[],params);
    prev_h = y(:,2);
    date = t(end);
    for i = 1:(length(t)-1)
       if abs((prev_h(i+1)-prev_h(i))/(t(i+1)-t(i))) < 1e-5
            if (t(i) > 5 || y(i,2) > 0.15 || y(i,2) < 1e-4)
                date = t(i);
                break
            end
       end
    end
    date = round(date);
end
function dx = SISV_fdem(t,x,z)
    dx = SISV_Model(t,x,[z 0.00004 0.1 0.00004 0.1]); 
end
function dx = SISV_C_fdem(t,x,z)
    dx = SISV_Model_Cumulative(t,x,[z 0.00004 0.1 0.00004 0.1]);
end