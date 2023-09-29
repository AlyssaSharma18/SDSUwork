model = @TuncerSIR_model;
params = [0.25 0.5];
y0 = [.99 .01 0];
error_levels = linspace(0,.3,31);
num_iters = 100;
init_guess = [0.1 0.1];
cumul = @cumulative;
tspan = linspace(0,50,51);
%% Changing Upper Bounds
ubs = [[Inf Inf];[2.5 5];[0.5 1];[3/8 3/4];1.1.*[0.25 0.5]];
ares = zeros(size(ubs,1),length(error_levels),2);
pes = zeros(size(ubs,1),length(error_levels),num_iters,2);
for i = 1:size(ubs,1)
    fmfun = fmsbnd_bnds([0 0],ubs(i,:));
    %Need to actually get fminfun
    [ares(i,:,:) pes(i,:,:,:)] = ARE(model,cumul,params,y0,init_guess,tspan,num_iters,error_levels,[],fmfun,[],[]);
    disp(i);
end
%% Plotting and Saving
% save('SIR_ARE_fmsb_bnds.mat','ares','pes');
error_levels_p = 100*error_levels;
figure
hold on
for i = 1:length(ubs)
    plot(error_levels_p,ares(i,:,1),'LineWidth',1.5)
end
hold off
legend('No Bounds','10x','2x','1.5x','1.1x');
xlabel('Error Level (%)'); ylabel('Average Relative Error'); title('ARE of \alpha for Upper Bounds of fminsearchbnd');
% saveas(gcf,'SIR_ARE_fmsb_a.png');
figure
hold on
for i = 1:length(ubs)
    plot(error_levels_p,ares(i,:,1),'LineWidth',1.5)
end
hold off
legend('No Bounds','10x','2x','1.5x','1.1x');
xlabel('Error Level (%)'); ylabel('Average Relative Error'); title('ARE of \beta for Upper Bounds of fminsearchbnd');
% saveas(gcf,'SIR_ARE_fmsb_b.png');
%% Fminfuns
function fmconbfun = fmcon_bnds(lb,ub)
    options = optimoptions('fmincon','Display','none');
    fmconbfun = @(error_fun,init_guess)(fmincon(error_fun,init_guess,[],[],[],[],lb,ub,[],options));
end
function fmsbndfun = fmsbnd_bnds(lb,ub)
    options = optimset('Display','none');
    fmsbndfun = @(error_fun,init_guess)(fminsearchbnd(error_fun,init_guess,lb,ub,options));
end
%% Function
function data = cumulative(y)
    start = y(1,2) + y(1,3);
    data = y(:,2)+y(:,3)-start;
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