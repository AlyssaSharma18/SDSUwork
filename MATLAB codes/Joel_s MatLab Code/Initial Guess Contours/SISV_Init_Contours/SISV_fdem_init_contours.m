% Will need to change for greater precision
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
init_guess = params; %Always have correct initial guess
y0 = [.99 .01 .99 0.01];
num_iters = 10;
error_levels = 0.1;
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
%% FMinSearch
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],[],[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Nelder-Mead');
% save('SISV_Init_bhbv_nm.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_nm.png');
% saveas(f2,'SISV_Init_bhbv_bv_nm.png');
% saveas(f3,'SISV_Init_bhbv_g_nm.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],[],[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Nelder-Mead');
% save('SISV_Init_bhg_nm.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_nm.png');
% saveas(f2,'SISV_Init_bhg_bv_nm.png');
% saveas(f3,'SISV_Init_bhg_g_nm.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],[],[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Nelder-Mead');
% save('SISV_Init_bvg_nm.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_nm.png');
% saveas(f2,'SISV_Init_bvg_bv_nm.png');
% saveas(f3,'SISV_Init_bvg_g_nm.png');
%% FMinCon
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmcon,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Interior-Point');
% save('SISV_Init_bhbv_ip.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_ip.png');
% saveas(f2,'SISV_Init_bhbv_bv_ip.png');
% saveas(f3,'SISV_Init_bhbv_g_ip.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmcon,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Interior-Point');
% save('SISV_Init_bhg_ip.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_ip.png');
% saveas(f2,'SISV_Init_bhg_bv_ip.png');
% saveas(f3,'SISV_Init_bhg_g_ip.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,2,48)];
gammas = [0 0.005 linspace(0.01,2,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmcon,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Interior-Point');
% save('SISV_Init_bvg_ip.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_ip.png');
% saveas(f2,'SISV_Init_bvg_bv_ip.png');
% saveas(f3,'SISV_Init_bvg_g_ip.png');
%% FMinUnc
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0,[],@fmunc,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Quasi-Newton');
% save('SISV_Init_bhbv_qn.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_qn.png');
% saveas(f2,'SISV_Init_bhbv_bv_qn.png');
% saveas(f3,'SISV_Init_bhbv_g_qn.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmunc,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Quasi-Newton');
% save('SISV_Init_bhg_qn.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_qn.png');
% saveas(f2,'SISV_Init_bhg_bv_qn.png');
% saveas(f3,'SISV_Init_bhg_g_qn.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmunc,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Quasi-Newton');
% save('SISV_Init_bvg_qn.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_qn.png');
% saveas(f2,'SISV_Init_bvg_bv_qn.png');
% saveas(f3,'SISV_Init_bvg_g_qn.png');
%% SQP
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmsqp,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','SQP');
% save('SISV_Init_bhbv_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_sqp.png');
% saveas(f2,'SISV_Init_bhbv_bv_sqp.png');
% saveas(f3,'SISV_Init_bhbv_g_sqp.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmsqp,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','SQP');
% save('SISV_Init_bhg_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_sqp.png');
% saveas(f2,'SISV_Init_bhg_bv_sqp.png');
% saveas(f3,'SISV_Init_bhg_g_sqp.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmsqp,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','SQP');
% save('SISV_Init_bvg_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_sqp.png');
% saveas(f2,'SISV_Init_bvg_bv_sqp.png');
% saveas(f3,'SISV_Init_bvg_g_sqp.png');
%% Active Set
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmactset,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Active Set');
% save('SISV_Init_bhbv_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_acts.png');
% saveas(f2,'SISV_Init_bhbv_bv_acts.png');
% saveas(f3,'SISV_Init_bhbv_g_acts.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmactset,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Active Set');
% save('SISV_Init_bhg_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_acts.png');
% saveas(f2,'SISV_Init_bhg_bv_acts.png');
% saveas(f3,'SISV_Init_bhg_g_acts.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmactset,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Active Set');
% save('SISV_Init_bvg_sqp.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_acts.png');
% saveas(f2,'SISV_Init_bvg_bv_acts.png');
% saveas(f3,'SISV_Init_bvg_g_acts.png');
%% Bounded interior point
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmcontenx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Bounded Interior Point');
% save('SISV_Init_bhbv_bip.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_btip.png');
% saveas(f2,'SISV_Init_bhbv_bv_btip.png');
% saveas(f3,'SISV_Init_bhbv_g_btip.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmactset,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Bounded Interior Point');
% save('SISV_Init_bhg_bip.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_btip.png');
% saveas(f2,'SISV_Init_bhg_bv_btip.png');
% saveas(f3,'SISV_Init_bhg_g_btip.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmactset,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Bounded Interior Point');
% save('SISV_Init_bvg_btip.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_btip.png');
% saveas(f2,'SISV_Init_bvg_bv_btip.png');
% saveas(f3,'SISV_Init_bvg_g_btip.png');
%% fminsearchbnd-bounded NM
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmsbndtenx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Bounded Nelder-Mead');
% save('SISV_Init_bhbv_btnm.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_btnm.png');
% saveas(f2,'SISV_Init_bhbv_bv_btnm.png');
% saveas(f3,'SISV_Init_bhbv_g_btnm.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmsbndtenx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Bounded Nelder-Mead');
% save('SISV_Init_bhg_btnm.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_btnm.png');
% saveas(f2,'SISV_Init_bhg_bv_btnm.png');
% saveas(f3,'SISV_Init_bhg_g_btnm.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmsbndtenx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Bounded Nelder-Mead');
% save('SISV_Init_bvg_btnm.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_btnm.png');
% saveas(f2,'SISV_Init_bvg_bv_btnm.png');
% saveas(f3,'SISV_Init_bvg_g_btnm.png');
%% fminsearchbnd-bounded NM at hundredx
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmsbndhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Bounded Nelder-Mead');
% save('SISV_Init_bhbv_bhnm.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_bhnm.png');
% saveas(f2,'SISV_Init_bhbv_bv_bhnm.png');
% saveas(f3,'SISV_Init_bhbv_g_bhnm.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmsbndhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Bounded Nelder-Mead');
% save('SISV_Init_bhg_bhnm.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_bhnm.png');
% saveas(f2,'SISV_Init_bhg_bv_bhnm.png');
% saveas(f3,'SISV_Init_bhg_g_bhnm.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmsbndhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Bounded Nelder-Mead');
% save('SISV_Init_bvg_bhnm.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_bhnm.png');
% saveas(f2,'SISV_Init_bvg_bv_bhnm.png');
% saveas(f3,'SISV_Init_bvg_g_bhnm.png');
%% fmincon-bounded IP at hundredx
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Bounded Interior Point');
% save('SISV_Init_bhbv_bhip.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_bhip.png');
% saveas(f2,'SISV_Init_bhbv_bv_bhip.png');
% saveas(f3,'SISV_Init_bhbv_g_bhip.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Bounded Interior Point');
% save('SISV_Init_bhg_bhip.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_bhip.png');
% saveas(f2,'SISV_Init_bhg_bv_bhip.png');
% saveas(f3,'SISV_Init_bhg_g_bhip.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Unbounded Nelder-Mead');
% save('SISV_Init_bvg_bhip.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_bip.png');
% saveas(f2,'SISV_Init_bvg_bv_bhip.png');
% saveas(f3,'SISV_Init_bvg_g_bhip.png');
%% fminsearchbnd no bounds
%% betah vs betav
beta_hs = [0 0.005 linspace(0.01,1,48)];
beta_vs = [0 0.005 linspace(0.01,1,48)];
[H,V] = meshgrid(beta_hs,beta_vs);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(beta_vs)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), V(r,c), gamma],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,V,[0.1 0.2],'\beta_h','\beta_v','Unbounded Nelder-Mead');
% save('SISV_Init_bhbv_bnnm.mat',"ares");
% saveas(f1,'SISV_Init_bhbv_bh_bnnm.png');
% saveas(f2,'SISV_Init_bhbv_bv_bnnm.png');
% saveas(f3,'SISV_Init_bhbv_g_bnnm.png');
%% Beta_h vs gamma
beta_hs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[H,G] = meshgrid(beta_hs,gammas);
ares = zeros([size(H),3]);
for c = 1:length(beta_hs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[H(r,c), beta_v, G(r,c)],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,H,G,[0.1 0.1],'\beta_h','\gamma','Unbounded Nelder-Mead');
% save('SISV_Init_bhg_bhip.mat',"ares");
% saveas(f1,'SISV_Init_bhg_bh_bhip.png');
% saveas(f2,'SISV_Init_bhg_bv_bhip.png');
% saveas(f3,'SISV_Init_bhg_g_bhip.png');
%% Beta_v vs gamma
beta_vs = [0 0.005 linspace(0.01,1,48)];
gammas = [0 0.005 linspace(0.01,1,48)];
[V,G] = meshgrid(beta_vs,gammas);
ares = zeros([size(G),3]);
for c = 1:length(beta_vs)
    for r = 1:length(gammas)
        ares(r,c,:) = ARE(model,prev_b,params,y0,[beta_h, V(r,c), G(r,c)],tspan,num_iters,0.1,[],@fmconhunx,[],@under_hundred);
    end
    fprintf('Column %i done\n',c);
end
%% Plotting
[f1,f2,f3] = Plot_Init_Contours(ares,V,G,[0.2 0.1],'\beta_v','\gamma','Unbounded Nelder-Mead');
% save('SISV_Init_bvg_bnnm.mat',"ares");
% saveas(f1,'SISV_Init_bvg_bh_bnnm.png');
% saveas(f2,'SISV_Init_bvg_bv_bnnm.png');
% saveas(f3,'SISV_Init_bvg_g_bnnm.png');
%% Functions
function [fig1,fig2,fig3] = Plot_Init_Contours(ares,X,Y,true_vals,xlab,ylab,algo_name)
    fig1 = figure;
    hold on
    contourf(X,Y,ares(:,:,1),[0 5 10 15 20 50]);
    scatter(true_vals(1),true_vals(2),250,'r','filled','pentagram');
    hold off
    xlabel(xlab);
    ylabel(ylab);
    title(['ARE of \beta_h using ', algo_name]);
    colorbar;
    fig2 = figure;
    hold on
    contourf(X,Y,ares(:,:,2),[0 5 10 15 20 50]);
    scatter(true_vals(1),true_vals(2),250,'r','filled','pentagram');
    hold off
    xlabel(xlab);
    ylabel(ylab);
    title(['ARE of \beta_v using ',algo_name]);
    colorbar;
    fig3 = figure;
    hold on
    contourf(X,Y,ares(:,:,3),[0 5 10 15 20 50]);
    scatter(true_vals(1),true_vals(2),250,'r','filled','pentagram');
    hold off
    xlabel(xlab);
    ylabel(ylab);
    title(['ARE of \gamma using ',algo_name]);
    colorbar;
end
function data = prev_both(y)
data = y(:,[2 4]);
end
%% Optimization Functions
function newpars = fmcon(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);
end
function newpars = fmunc(error_fun,init_guess)
    options = optimoptions('fminunc','Display','none');
    newpars = fminunc(error_fun,init_guess,options);
end
function newpars = fmsqp(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none','Algorithm','sqp');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);
end
function newpars = fmactset(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none','Algorithm','active-set');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);
end
function newpars = fmcontenx(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[0 0 0],[1 2 1],[],options);
end
function newpars = fmsbndtenx(error_fun,init_guess)
    options = optimset('Display','none');
    newpars = fminsearchbnd(error_fun,init_guess,[0 0 0],[1 2 1],options);
end
function newpars = fmconhunx(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[0 0 0],[10 20 10],[],options);
end
function newpars = fmsbndhunx(error_fun,init_guess)
    options = optimset('Display','none');
    newpars = fminsearchbnd(error_fun,init_guess,[0 0 0],[10 20 10],options);
end
function newpars = fmsbnd(error_fun,init_guess)
    options = optimset('Display','none');
    newpars = fminsearchbnd(error_fun,init_guess,[0 0 0],[Inf Inf Inf],options);
end
function date = end_date(model,y0,params)
    fulltspan = linspace(0,10000,10000);
    [t,y] = ode45(model,fulltspan,y0,[],params);
    prev_h = y(:,2);
    date = t(end);
    for i = 1:(length(t)-1)
       if ((prev_h(i+1)-prev_h(i))/(t(i+1)-t(i)) < 1e-7)
           date = t(i);
           break
       end
    end
    date = round(date);
end
function is = isStiff(params)
    r0 = sqrt(params(1)*params(2)/((0.00004 + params(3))*0.1));
    if (r0 > 50)
        is = true;
    else
        is = false;
    end
end
function is = under_hundred(params)
    if (params(1) > 100) || (params(2) > 100)
        is = true;
    else
        is = false;
    end
end
function dx = SISV_fdem(t,x,z)
    dx = SISV_Model(t,x,[z 0.00004 0.1 0.00004 0.1]); 
end