model = @fixedpimu;
prev = @prev_both;
params = [0.0001 0.001 0.08996];
y0 = [99 1 999 1];
num_iters = 5;
error_levels = 0;
tspan = linspace(0,1000,101);
%% Running for fminsearch
[betah_g,betav_g] = meshgrid(0.00001:0.00001:0.0002,0.0001:0.0001:0.002);
s = size(betah_g);
cols = s(2);
rows = s(1);
ares = zeros([s 3]);
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[betah_g(r,c) betav_g(r,c) 0.08996],tspan,num_iters,error_levels,[],[]);
    end
end
%% Plotting
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
title('ARE of \beta_v for Initial Conditions at \sigma=0 using fminsearch');
xlabel('\beta_h');
ylabel('\beta_v');
writematrix(ares,'SIV_AREs_fdem_Fminsearch_contour.xlsx');
colorbar;
% saveas(gcf,'SIV_init_fdem_fmsearch_betah.png');
colorbar;
%% Running for fmincon
[betah_g,betav_g] = meshgrid(0.00001:0.00001:0.0002,0.0001:0.0001:0.002);
s = size(betah_g);
cols = s(2);
rows = s(1);
ares = zeros([s 3]);
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[betah_g(r,c) betav_g(r,c) 0.08996],tspan,num_iters,error_levels,[],@fmcon);
    end
end
%% Plotting
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
title('ARE of \beta_v for Initial Conditions at \sigma=0 using fmincon');
xlabel('\beta_h');
ylabel('\beta_v');
writematrix(ares,'SIV_AREs_fdem_Fmincon_contour.xlsx');
colorbar;
% saveas(gcf,'SIV_init_fdem_fmcon_betah.png');
%% Running for fminunc
[betah_g,betav_g] = meshgrid(0.00001:0.00001:0.0002,0.0001:0.0001:0.002);
s = size(betah_g);
cols = s(2);
rows = s(1);
ares = zeros([s 3]);
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[betah_g(r,c) betav_g(r,c) 0.08996],tspan,num_iters,error_levels,[],@fmunc);
    end
end
%% Plotting
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
title('ARE of \beta_v for Initial Conditions at \sigma=0 using fminunc');
xlabel('\beta_h');
ylabel('\beta_v');
writematrix(ares,'SIV_AREs_fdem_Fminunc_contour.xlsx');
colorbar;
% saveas(gcf,'SIV_init_fdem_fmunc_betah.png');
%% Functions
function data = prev_both(y)
data = y(:,[2 4]);
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