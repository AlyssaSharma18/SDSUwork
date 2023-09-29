model = @SISV_Model;
prev_b = @prev_both;
pih = 0.00004;
piv = 0.1;
muh = 0.00004;
muv = 0.1;
betah = 0.1;
betav = 0.2;
gamma = .09996;
params = [betah betav gamma muh muv pih piv];
tspan = linspace(0,216,216);
y0 = [.99 .01 .99 .01];
[~,true_y] = ode45(model,tspan,y0,[],params);
true_data = prev_b(true_y);
%% Running betah vs betav
beta_hs = linspace(0.01,0.2,250);
beta_vs = linspace(0.01,0.4,250);
[H,V] = meshgrid(beta_hs,beta_vs);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) V(r,c) gamma muh muv pih piv],y0,tspan,@OLS);
    end
end
%% Plotting
figure
hold on
contourf(H,V,errors,[0 0.25 0.75 1.5 5 10 19 100 1e3]);
plot(beta_hs,(0.02)./beta_hs,'m','LineWidth',1.5);
scatter(betah,betav,250,'r','filled','pentagram');
ylim([0.01 0.4]);
hold off
xlabel('\beta_h');
ylabel('\beta_v');
title('Error For \beta_h vs \beta_v Using OLS')
colorbar;
set(gca,'ColorScale','log');
%% Running betav vs gamma
beta_vs = linspace(0.01,0.4,250);
gammas = linspace(0.01,0.2,250);
[V,G] = meshgrid(beta_vs,gammas);
s = size(V);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[betah V(r,c) G(r,c) muh muv pih piv],y0,tspan,@OLS);
    end
end
%% Plotting
figure
hold on
contourf(V,G,errors,[0 0.5 1 2.5 5 10 19 100 1000]);
plot(beta_vs,(beta_vs-(2*muh))./2,'m','LineWidth',1.5);
scatter(betav,gamma,250,'r','filled','pentagram');
xlabel('\beta_v');
ylabel('\gamma');
ylim([0.01,0.2]);
title('Error For \beta_v vs \gamma Using OLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betah vs gamma
beta_hs = linspace(0.01,0.2,250);
gammas = linspace(0.01,0.2,250);
[H,G] = meshgrid(beta_hs,gammas);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) betav G(r,c) muh muv pih piv],y0,tspan,@OLS);
    end
end
%% Plotting
figure
hold on
contourf(H,G,errors,[0 0.1 0.5 2 5 10 19 100 1e4]);
plot(beta_hs,(beta_hs-(2*muh)),'m','LineWidth',1.5);
scatter(betah,gamma,250,'r','filled','pentagram');
hold off
xlabel('\beta_h')
ylabel('\gamma');
ylim([0.01,0.2]);
title('Error For \beta_h vs \gamma Using OLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betav vs muv
beta_vs = linspace(0.01,0.4,250);
mu_vs = linspace(0.01,0.2,250);
[V,M] = meshgrid(beta_vs,mu_vs);
s = size(V);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[betah V(r,c) gamma muh M(r,c) pih M(r,c)],y0,tspan,@OLS);
    end
end
%% Plotting
figure
hold on
contourf(V,M,errors,[0,0.1,1,2,5,10,19,100]);
plot(beta_vs,(beta_vs./2),'m','LineWidth',1.5);
scatter(betav,muv,250,'r','filled','pentagram');
ylim([0.01 0.2]);
hold off
xlabel('\beta_v');
ylabel('\mu_v');
title('Error For \beta_v vs \mu_v Using OLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betah vs muh
beta_hs = linspace(0.01,0.2,1000);
mu_hs = linspace(0,0.00008,100);
[H,M] = meshgrid(beta_hs,mu_hs);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) betav gamma M(r,c) muv M(r,c) piv],y0,tspan,@OLS);
    end
end
%% Plotting
% Log scale somewhat misses the pointyness
figure
hold on
contourf(H,M,errors,[0,0.1,1,5,10,19,100]);
plot(beta_hs,beta_hs-gamma,'m','LineWidth',1.5);
scatter(betah,muh,250,'r','filled','pentagram');
hold off
xlabel('\beta_h');
ylabel('\mu_h');
ylim([0 0.00008]);
title('Error For \beta_h vs \mu_h Using OLS');
colorbar;
set(gca,'ColorScale','log');
%% Now All with GLS





%% Running betah vs betav
beta_hs = linspace(0,0.2,250);
beta_vs = linspace(0,0.4,250);
[H,V] = meshgrid(beta_hs,beta_vs);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) V(r,c) gamma muh muv pih piv],y0,tspan,@GLS);
    end
end
%% Plotting
figure
hold on
contourf(H,V,errors,[0 10 25 50 100 400 1000 10000]);
plot(beta_hs,(0.02)./beta_hs,'r','LineWidth',1.5);
scatter(betah,betav,250,'r','filled','pentagram');
ylim([0 0.4]);
hold off
xlabel('\beta_h');
ylabel('\beta_v');
title('Error For \beta_h vs \beta_v Using GLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betav vs gamma
beta_vs = linspace(0,0.4,100);
gammas = linspace(0,0.2,100);
[V,G] = meshgrid(beta_vs,gammas);
s = size(V);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[betah V(r,c) G(r,c) muh muv pih piv],y0,tspan,@GLS);
    end
end
%% Plotting
figure
hold on
contourf(V,G,errors,[0 5 20 50 100 1000 1e4]);
plot(beta_vs,(beta_vs-(2*muh))./2,'m','LineWidth',1.5);
scatter(betav,gamma,250,'r','filled','pentagram');
xlabel('\beta_v');
ylabel('\gamma');
ylim([0,0.2]);
title('Error For \beta_v vs \gamma Using GLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betah vs gamma
beta_hs = linspace(0.01,0.2,250);
gammas = linspace(0.01,0.2,250);
[H,G] = meshgrid(beta_hs,gammas);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) betav G(r,c) muh muv pih piv],y0,tspan,@GLS);
    end
end
%% Plotting
figure
hold on
contourf(H,G,errors,[0 2.5 25 100 1000 1e4]);
plot(beta_hs,(beta_hs-(2*muh)),'m','LineWidth',1.5);
scatter(betah,gamma,250,'r','filled','pentagram');
hold off
xlabel('\beta_h')
ylabel('\gamma');
title('Error For \beta_h vs \gamma Using GLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betav vs muv
beta_vs = linspace(0.01,0.4,250);
mu_vs = linspace(0.01,0.2,250);
[V,M] = meshgrid(beta_vs,mu_vs);
s = size(V);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[betah V(r,c) gamma muh M(r,c) pih M(r,c)],y0,tspan,@GLS);
    end
end
%% Plotting
figure
hold on
contourf(V,M,errors,[0 2.5 25 50 100 1000 1e4]);
plot(beta_vs,(beta_vs./2),'m','LineWidth',1.5);
scatter(betav,muv,250,'r','filled','pentagram');
ylim([0.01 0.2]);
hold off
xlabel('\beta_v');
ylabel('\mu_v');
title('Error For \beta_v vs \mu_v Using GLS');
colorbar;
set(gca,'ColorScale','log');
%% Running betah vs muh
beta_hs = linspace(0.01,0.2,1000);
mu_hs = linspace(0,0.00008,100);
[H,M] = meshgrid(beta_hs,mu_hs);
s = size(H);
cols = s(2);
rows = s(1);
errors = zeros(s);
parfor c = 1:cols
    for r = 1:rows
        errors(r,c) = Model_Error(model,prev_b,true_data,[H(r,c) betav gamma M(r,c) muv M(r,c) piv],y0,tspan,@GLS);
    end
end
%% Plotting
figure
hold on
contourf(H,M,errors,[0,0.05,5,100,250,500,1000,1e4]);
plot(beta_hs,beta_hs-gamma,'m','LineWidth',1.5);
scatter(betah,muh,250,'r','filled','pentagram');
hold off
xlabel('\beta_h');
ylabel('\mu_h');
ylim([0,0.00008]);
title('Error For \beta_h vs \mu_h Using GLS');
colorbar;
set(gca,'ColorScale','log');
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
function err = Model_Error(model,observations,true_data,params,y0,tspan,obj_fun)
    %Checks for negative parameters
    if any(params<0)
        err = 100000; %Try inf?
    else
        [~,y_new] = ode45(model,tspan,y0,[],params);
        data_new = observations(y_new);
        err = obj_fun(true_data,data_new);
    end
end
function er = OLS(data1,data2)
    diffsql2 = sum((data1 - data2).^2,2);
    er = sum(diffsql2);%Use dimensions for preciseness instead of this nonsense (l2 norm)
end
function er = GLS(data1,data2)
    reldiffsq = ((data1-data2).^2)./(data1.^2);
    er = sum(reldiffsq,"all");
end