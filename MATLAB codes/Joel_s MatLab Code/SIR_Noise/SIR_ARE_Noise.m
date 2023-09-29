%% Setup
error_levels = linspace(0,.3,10);
init_guess = [0.1 0.1];
model = @TuncerSIR_model;
prev = @prevalence;
params = [0.25 0.5];
y0 = 100.*[0.99 0.01 0];
num_iters = 25;
tspan = linspace(0,50,51);
[~,y] = ode45(model,tspan,y0,[],params);
figure
plot(tspan,y(:,2))
%% For Relative Noise
noisemaker = @(sigma,data)(Weight_Error(sigma,data,1,1));
ols = @(data1,data2,true_data)(WLS(data1,data2,1,0,true_data));
% poisson = @(data1,data2,true_data)(WLS(data1,data2,1,1/2,true_data));
gls = @(data1,data2,true_data)(WLS(data1,data2,1,1,true_data));
obj_funs = {ols,gls};
ares = zeros([length(obj_funs),length(error_levels),2]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker);
end
%% Plotting
error_levels_p = 100.*error_levels;
% save('SIR_ARE_Rel_Noise.mat','ares');
%beta
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
% plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
legend('OLS','GLS');ylim([0 min(yl(2),100)]);
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \beta for Relative Noise');
% saveas(gcf,'SIR_ARE_Rel_Noise_a.png')
%alpha
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
% plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \alpha for Relative Noise');
legend('OLS','GLS');ylim([0 min(yl(2),100)]);
% saveas(gcf,'SIR_ARE_Rel_Noise_b.png')
%% For Poisson Noise
noisemaker = @(sigma,data)(Weight_Error(sigma,data,1,0.5));
ols = @(data1,data2,true_data)(WLS(data1,data2,1,0,true_data));
poisson = @(data1,data2,true_data)(WLS(data1,data2,1,1/2,true_data));
gls = @(data1,data2,true_data)(WLS(data1,data2,1,1,true_data));
obj_funs = {ols,poisson,gls};
ares = zeros([3,length(error_levels),2]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker);
end
%% Plotting
error_levels_p = 100.*error_levels;
%beta
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \beta for Poisson Noise');
legend('OLS','Poisson','GLS')
ylim([0 min(yl(2),100)]);
%alpha
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \alpha for Poisson Noise');
legend('OLS','Poisson','GLS');
ylim([0 min(yl(2),100)]);
%% For Absolute Noise (10% of pop)
noisemaker = @(sigma,data)(Weight_Error(sigma,data,100,0));
ols = @(data1,data2,true_data)(WLS(data1,data2,1,0,true_data));
% poisson = @(data1,data2,true_data)(WLS(data1,data2,1,1/2,true_data));
gls = @(data1,data2,true_data)(WLS(data1,data2,1,1,true_data));
obj_funs = {ols,gls};
ares = zeros([length(obj_funs),length(error_levels),2]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker);
end
%% Plotting Absolute Noise
error_levels_p = 100.*error_levels;
save('SIR_ARE_Abs_Noise.mat','ares');
%beta
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
% plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
legend('OLS','GLS');ylim([0 min(yl(2),100)]);
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \beta for Absolute Noise');
saveas(gcf,'SIR_ARE_Abs_Noise_a.png')
%alpha
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
% plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)'); ylabel('Average Relative Error');title('ARE of \alpha for Absolute Noise');
legend('OLS','GLS');ylim([0 min(yl(2),100)]);
saveas(gcf,'SIR_ARE_Abs_Noise_b.png')
%% Weighted before Peak
weight = linspace(0,50,51);
error_levels = linspace(0,.3,50);
noisemaker = @(sigma,data)(Weight_Error(sigma,data,1,1));
for i = 1:length(weight)
    if tspan(i) < 17
        weight(i) = 2;
    else
        weight(i) = 0.5;
    end
end
ols = @(data1,data2)(WLS(data1,data2,1,0));
poisson = @(data1,data2)(WLS(data1,data2,1,1/2));
gls = @(data1,data2)(WLS(data1,data2,1,1));
ols_w = @(data1,data2)(WLS(data1,data2,weight,0));
poisson_w = @(data1,data2)(WLS(data1,data2,weight,1/2));
gls_w = @(data1,data2)(WLS(data1,data2,weight,1));
obj_funs = {ols,ols_w,poisson,poisson_w,gls,gls_w};
ares = zeros([6,length(error_levels),2]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker);
end
%% Plotting
figure
hold on
plot(error_levels,ares(1,:,1),'b','LineWidth',1.5);
plot(error_levels,ares(2,:,1),'b--','LineWidth',1.5);
plot(error_levels,ares(3,:,1),'r','LineWidth',1.5);
plot(error_levels,ares(4,:,1),'r--','LineWidth',1.5);
plot(error_levels,ares(5,:,1),'m','LineWidth',1.5);
plot(error_levels,ares(6,:,1),'m--','LineWidth',1.5);
hold off
legend('OLS','OLS_W','Poisson','Poisson_W','GLS','GLS_W');
xlabel('Error Levels');
ylabel('Average Relative Error');
title('ARE of \alpha for Objective Functions using Weighted Data')
figure
hold on
plot(error_levels,ares(1,:,2),'b','LineWidth',1.5);
plot(error_levels,ares(2,:,2),'b--','LineWidth',1.5);
plot(error_levels,ares(3,:,2),'r','LineWidth',1.5);
plot(error_levels,ares(4,:,2),'r--','LineWidth',1.5);
plot(error_levels,ares(5,:,2),'m','LineWidth',1.5);
plot(error_levels,ares(6,:,2),'m--','LineWidth',1.5);
legend('OLS','OLS_W','Poisson','Poisson_W','GLS','GLS_W');
xlabel('Error Levels');
ylabel('Average Relative Error');
title('ARE of \beta for Objective Functions using Weighted Data')
hold off
%% Old
% %% Basic ARE (for comparison)
% are_noweight = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels);
% %% Around 10 (test)
% weight = ones(length(tspan),1);
% for i = 1:length(weight)
%     if abs(tspan(i)-10) < 17
%         weight(i) = 2;
%     else
%         weight(i) = 0.5;
%     end
% end
% are_weight = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[],@(sigma,data)(Weight_Relative(sigma,data,weight)));
% are_adjusted = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,@(data1,data2)(WLS(data1,data2,weight)),[],@(sigma,data)(Weight_Relative(sigma,data,weight)));
% %% Plot
% figure
% hold on
% plot(error_levels,are_noweight(:,1),'b','LineWidth',1.5);
% plot(error_levels,are_noweight(:,2),'b--','LineWidth',1.5);
% plot(error_levels,are_weight(:,1),'r','LineWidth',1.5);
% plot(error_levels,are_weight(:,2),'r--','LineWidth',1.5);
% plot(error_levels,are_adjusted(:,1),'m','LineWidth',1.5);
% plot(error_levels,are_adjusted(:,2),'m--','LineWidth',1.5);
% hold off
% legend('No Weight:\alpha','No Weight:\beta','Weight:\alpha','Weight:\beta');
%% Functions
function noisy_data = Weight_Error(sigma,data,weight,p)
    noisy_data = data + weight.*(data.^p).*normrnd(0,sigma,size(data));
end
function fval = WLS(data1,data2,weight,p,true_data)
    if isempty(true_data)
        true_data = data1;
        disp('Using Noisy Data for Weighting');
    end
    vector = (weight.^(-2))'.*(abs(true_data).^(-2*p)).*((data1-data2).^2);
    fval = sum(vector);
end
function data = prevalence(y)
data = y(:,2);
end