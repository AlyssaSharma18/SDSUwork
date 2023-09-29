% Cumulative Data for both, with cumulative error? With different bounds for fmincon,etc.
%% Setup
model = @SISV_fdem_Cumulative;
betah = 0.1;
betav = 0.2;
gamma = 0.1;
params = [betah betav gamma];
y0 = [0.99 0.01 0.99 0.01 0 0];
e_date = end_date(model,y0,params);
tspan = linspace(0,216,217);
cumul_b = @(y)(y(:,[5 6]));
error_levels = linspace(0,.3,10);
num_iters = 10;
%% Changing fmsbnd's bounds
%Init guess is 0.05s
init_guess = [0.05 0.1 0.05];
prev_b = @(y)(y([2 4]));
ubs = [[Inf Inf Inf];[1 2 1];[0.5 1 0.5];[0.2 0.4 0.2];[0.15 0.3 0.15];[0.11 0.22 0.11]];
ares = zeros(size(ubs,1),length(error_levels),3);
pes = zeros(size(ubs,1),length(error_levels),num_iters,3);
for i = 1:size(ubs,1)
    fmfun = fmsbnd_bnds([0 0 0],ubs(i,:));
    [ares(i,:,:) pes(i,:,:,:)] = ARE(model,cumul_b,params,y0,init_guess,tspan,num_iters,error_levels,[],fmfun,[],@under_hundred);
    disp(i);
end
%% Plotting and saving
% save('SISV_ARE_fmsb_bnds.mat','ares','pes');
error_levels_p = 100*error_levels;
figure
hold on
for i = 1:length(ubs)
    plot(error_levels_p,ares(i,:,1),'LineWidth',1.5)
end
hold off
legend('No Bounds','10x','5x','2x','1.5x','1.1x');
xlabel('Error Level (%)'); ylabel('Average Relative Error'); title('ARE of \beta_h for Upper Bounds of fminsearchbnd');
% saveas(gcf,'SISV_ARE_fmsb_bh.png');
figure
hold on
for i = 1:length(ubs)
    plot(error_levels_p,ares(i,:,2),'LineWidth',1.5)
end
hold off
legend('No Bounds','10x','5x','2x','1.5x','1.1x');
xlabel('Error Level (%)'); ylabel('Average Relative Error'); title('ARE of \beta_v for Upper Bounds of fminsearchbnd');
% saveas(gcf,'SISV_ARE_fmsb_bh.png');
figure
hold on
for i = 1:length(ubs)
    plot(error_levels_p,ares(i,:,3),'LineWidth',1.5)
end
hold off
legend('No Bounds','10x','5x','2x','1.5x','1.1x');
xlabel('Error Level (%)'); ylabel('Average Relative Error'); title('ARE of \gamma for Upper Bounds of fminsearchbnd');
% saveas(gcf,'SISV_ARE_fmsb_bh.png');
%%
error_levels = [0 0.01 0.05 0.1 0.2 0.3];
ub2 = [1 1 1];
lb2 = [0 0 0];
ub3 = [0.5 0.5 0.5];
lb3 = [0.01 0.01 0.01];
ub4 = [0.2 0.4 0.2];
lb4 = [0.02 0.02 0.02];
bnd_funs = {fmsbnd_bnds([],[]),fmsbnd_bnds(lb4,ub4),fmsbnd_bnds([0.04 0.04 0.04],[0.15 0.3 0.15])};
for i = 1:length(bnd_funs)
    are = ARE(model,cumul_b,params,y0,[0.05 0.05 0.05],tspan,num_iters,error_levels,[],bnd_funs{i},[],@under_hundred);
    disp(i);
    disp(are);
end
%% Functions
function dx = SISV_fdem(t,x,z)
    dx = SISV_Model(t,x,[z 0.00004 0.1 0.00004 0.1]); 
end
function fmconbfun = fmcon_bnds(lb,ub)
    options = optimoptions('fmincon','Display','none');
    fmconbfun = @(error_fun,init_guess)(fmincon(error_fun,init_guess,[],[],[],[],lb,ub,[],options));
end
function fmsbndfun = fmsbnd_bnds(lb,ub)
    options = optimset('Display','none');
    fmsbndfun = @(error_fun,init_guess)(fminsearchbnd(error_fun,init_guess,lb,ub,options));
end
function dx = SISV_fdem_Cumulative(t,x,z)
    dx = SISV_Model_Cumulative(t,x,[z(1) z(2) z(3) 0.00004 0.1 0.00004 0.1]);
end
function date = end_date(model,y0,params)
    fulltspan = linspace(0,10000,50000);
    [t,y] = ode45(model,fulltspan,y0,[],params);
    prev_h = y(:,2);
    prev_v = y(:,4);
    date = t(end);
    for i = 1:(length(t)-1)
       if (abs((prev_h(i+1)-prev_h(i))/(t(i+1)-t(i))) < 1e-6)
            if (t(i) > 1)
                date = t(i);
                break
            end
       end
    end
    date = round(date);
end
function is = under_hundred(params)
    if (params(1) > 100 || params(2) > 100)
        is = true;
    else
        is = false;
    end
end