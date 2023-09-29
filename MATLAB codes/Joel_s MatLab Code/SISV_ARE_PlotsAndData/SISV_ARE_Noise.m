%% Setup
model = @SISV_fdem;
error_levels = linspace(0,.3,31);
init_guess = [0.05 0.1 0.05];
prev_b = @(y)(y(:,[2 4]));
params = [0.1 0.2 0.1];
y0 = [99 1 99900 100];
num_iters = 25;
tspan = linspace(0,216,217);
%%
[t,y] = ode45(model,tspan,y0,[],params);
figure
plot(t,(y(:,[2 4]))./[100 100000]);
%Vpop = 100,000, hpop = 100
%% GOOD RESULT ALERT
%% With GLS Noise
error_levels = linspace(0,0.3,10);
noisemaker = @(sigma,data)(Weight_Error(sigma,data,1,1));
ols = @(data1,data2,true_data)(WLS(data1,data2,1,0,true_data));
ols_pop_adj = @(data1,data2,true_data)(WLS(data1,data2,[1;1000],0,true_data));
gls = @(data1,data2,true_data)(WLS(data1,data2,1,1,true_data));
obj_funs = {ols,ols_pop_adj,gls};
ares = zeros([length(obj_funs),length(error_levels),3]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker,@under_hundred);
end
%% Plotting
save('SISV_ARE_Rel_Noise.mat','ares');
error_levels_p = 100.*error_levels;
%betah
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \beta_h for Relative Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
saveas(gcf,'SISV_ARE_Rel_Noise_bh.png');
%betav
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \beta_v for Relative Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
saveas(gcf,'SISV_ARE_Rel_Noise_bv.png');
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \gamma for Relative Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
saveas(gcf,'SISV_ARE_Rel_Noise_g.png');
%% With OLS Noise
noisemaker = @(sigma,data)(Weight_Error(sigma,data,[10 10000],0));
ols = @(data1,data2,true_data)(WLS(data1,data2,1,0,true_data));
ols_pop_adj = @(data1,data2,true_data)(WLS(data1,data2,[1;1000],0,true_data));
gls = @(data1,data2,true_data)(WLS(data1,data2,1,1,true_data));
obj_funs = {ols,ols_pop_adj,gls};
ares = zeros([length(obj_funs),length(error_levels),3]);
for i = 1:length(obj_funs)
    ares(i,:,:) = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,obj_funs{i},[],noisemaker,@under_hundred);
end
%% Plotting
error_levels_p = 100.*error_levels;
%betah
figure
hold on
plot(error_levels_p,ares(1,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,1),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,1),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \beta_h for Absolute Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
%betav
figure
hold on
plot(error_levels_p,ares(1,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,2),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,2),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \beta_v for Absolute Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
%gamma
figure
hold on
plot(error_levels_p,ares(1,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(2,:,3),'LineWidth',1.5);
plot(error_levels_p,ares(3,:,3),'LineWidth',1.5);
yl = ylim;
plot(error_levels_p,error_levels_p,'--');
hold off
xlabel('Error Level (%)');ylabel('Average Relative Error');
title('ARE of \gamma for Absolute Noise with Different Objective Functions');
legend('OLS','Population Adjusted OLS','GLS');
ylim([0 min(100,yl(2))]);
%% Testing
noisemaker = @(sigma,data)(Weight_Error(sigma,data,1,1));
ols = @(data1,data2)(WLS(data1,data2,1,0));
gls = @(data1,data2)(WLS(data1,data2,1,1));
ols_pop_adj = @(data1,data2,true_data)(WLS(data1,data2,[100;100000],0,true_data));
data = [0.1 1;2 3; 20 21; 200 400];
noisemaker(0.1,data);
data2 = [0.1 1;2 3;20 21;210 410];
data3 = [0.1 1;1 3;20 21;200 400];
ols_pop_adj(data,data2,data)
%% Noise-Related Functions
function noisy_data = Weight_Error(sigma,data,weight,p)
    noisy_data = data + weight.*(data.^p).*normrnd(0,sigma,size(data));
end
function fval = WLS(data1,data2,weight,p,true_data)
    if isempty(true_data)
        true_data = data1;
        disp('Using Noisy Data for Weighting');
    end
    vector = (weight.^(-2))'.*(abs(true_data).^(-2*p)).*((data1-data2).^2);
    fval = sum(vector,'all');
end
%% Functions
function is = under_hundred(params)
    if (params(1) > 100 || params(2) > 100)
        is = true;
    else
        is = false;
    end
end
function date = end_date(model,y0,params)
    fulltspan = linspace(0,10000,50000);
    [t,y] = ode45(model,fulltspan,y0,[],params);
    prev_h = y(:,2);
    prev_v = y(:,4);
    date = t(end);
    for i = 1:(length(t)-1)
       if (abs((prev_h(i+1)-prev_h(i))/(t(i+1)-t(i))) < 1e-5) && (abs((prev_v(i+1)-prev_v(i))/(t(i+1)-t(i))))
            if (t(i) > 1 || prev_h(i) > 0.15 || prev_h(i) < 1e-3)
                date = t(i);
                break
            end
       end
    end
    date = round(date);
end
function dx = SISV_fdem(t,x,z)
    dx = SISV_Model(t,x,[z 0.00004 0.1 0.004 10000]); 
end