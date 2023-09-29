%% R0 = sqrt(2) Setup Pop of 1
% Note that Model does not change
%Start with 1/100 infected for each
model = @SISV_Model;
betah = 0.1;
betav = 0.2;
gamma = 0.09996;
muh = 0.00004;
muv = 0.1;
pih = 0.00004;
piv = 0.1;
params = [betah betav gamma muh muv pih piv];
y0 = [.99 0.01 .99 0.01];
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
prev_b = @prev_both;
error_levels = [0 0.01 0.05 0.1 0.2 0.3];
%% Fixing Nothing (sqrt(2))
num_iters = 50;
init_guess = [0.1 0.2 0.09996 0.00004 0.1 0.00004 0.1];
[are_c,pe_c] = ARE(model,prev_b,params,y0,init_guess,tspan,num_iters,error_levels,[],[],[],@under_hundred);
disp(are_c);
init_guess2 = [0.05 0.05 0.05 0.00004 0.1 0.00004 0.1];
[are_d,pe_d] = ARE(model,prev_b,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
disp(are_d);
save('SISV_fnone_AREs.mat','are_c','are_d','pe_c','pe_d');
matrix2latex(are_c,'fnone_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_h$','$\mu_v$','$\Pi_h$','$\Pi_v$'},'alignment','c','format','%.2f');
matrix2latex(are_d,'fnone_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_h$','$\mu_v$','$\Pi_h$','$\Pi_v$'},'alignment','c','format','%.2f');
%% Test
fig1 = Plot_Param_Estimates([pe_a(4,:,1);pe_a(4,:,2)],model,prev_b,[0.1 0.2],y0,tspan,@(p)([p(1) p(2) 0.1 0.00004 0.1 0.00004 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\beta_v');
title('Parameter Estimates Guess')
fig2 = Plot_Param_Estimates([pe_a(4,:,1);pe_a(4,:,3)],model,prev_b,[0.1 0.1],y0,tspan,@(p)([p(1) 0.2 p(2) 0.00004 0.1 0.00004 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\gamma');
title('Parameter Estimates Guess')
fig3 = Plot_Param_Estimates([pe_a(4,:,2);pe_a(4,:,3)],model,prev_b,[0.2 0.1],y0,tspan,@(p)([0.1 p(1) p(2) 0.00004 0.1 0.00004 0.1]),@OLS);
xlabel('\beta_v');
ylabel('\gamma');
title('Parameter Estimates Guess')
%% Fixing Populations
model = @SISV_fpops;
num_iters = 50;
params = [betah betav gamma muh muv];
init_guess1 = params;
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
init_guess2 = [0.05 0.05 0.05 0.00004 0.05];
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
save('SISV_fpops_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fpops_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_h$','$\mu_v$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fpops_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_h$','$\mu_v$'},'alignment','c','format','%.2f');
%% Plotting
fig1 = Plot_Param_Estimates([pe2(4,:,1);pe2(4,:,2)],model,prev_b,[0.1,0.2],y0,tspan,@(p)([p(1) p(2) 0.1 0.00004 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\beta_v');
% saveas(fig1,'SISV_PE_bhbv_guess_off.png');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
fig2 = Plot_Param_Estimates([pe2(4,:,1);pe2(4,:,3)],model,prev_b,[0.1,0.1],y0,tspan,@(p)([p(1) 0.2 p(2) 0.00004 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig2,'SISV_PE_bhg_guess_off.png');
fig3 = Plot_Param_Estimates([pe2(4,:,2);pe2(4,:,3)],model,prev_b,[0.2,0.1],y0,tspan,@(p)([0.1 p(1) p(2) 0.00004 0.1]),@OLS);
xlabel('\beta_v');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig3,'SISV_PE_bvg_guess_off.png');
fig4 = Plot_Param_Estimates([pe2(4,:,2);pe2(4,:,5)],model,prev_b,[0.2,0.1],y0,tspan,@(p)([0.1 p(1) 0.1 0.00004 p(2)]),@OLS);
xlabel('\beta_v');
ylabel('\mu_v');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig4,'SISV_PE_bvm_guess_off.png');
%% Fixing Host Demographics, Vector Population
model = @SISV_fh;
num_iters = 10;
params = [betah betav gamma muv];
init_guess1 = params;
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
init_guess2 = [0.05 0.05 0.05 0.05];
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
% save('SISV_fhdemvpop_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fhdemvpop_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_v$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fhdemvpop_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$','$\mu_v$'},'alignment','c','format','%.2f');
%% Plotting
fig1 = Plot_Param_Estimates([pe1(4,:,1);pe1(4,:,2)],model,prev_b,[0.1,0.2],y0,tspan,@(p)([p(1) p(2) 0.1 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\beta_v');
% saveas(fig1,'SISVfh_PE_bhbv_guess_off.png');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
fig2 = Plot_Param_Estimates([pe1(4,:,1);pe1(4,:,3)],model,prev_b,[0.1,0.1],y0,tspan,@(p)([p(1) 0.2 p(2) 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig2,'SISVfh_PE_bhg_guess_off.png');
fig3 = Plot_Param_Estimates([pe1(4,:,2);pe1(4,:,3)],model,prev_b,[0.2,0.1],y0,tspan,@(p)([0.1 p(1) p(2) 0.1]),@OLS);
xlabel('\beta_v');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig3,'SISVfh_PE_bvg_guess_off.png');
fig4 = Plot_Param_Estimates([pe1(4,:,2);pe1(4,:,4)],model,prev_b,[0.2,0.1],y0,tspan,@(p)([0.1 p(1) 0.1 p(2)]),@OLS);
xlabel('\beta_v');
ylabel('\mu_v');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig4,'SISVfh_PE_bvm_guess_off.png');
%% Fixing all Demographics
model = @SISV_fdem;
num_iters = 10;
params = [betah betav gamma];
init_guess1 = params;
init_guess2 = [0.05 0.05 0.05];
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
save('SISV_fdems_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fdems_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fdems_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
%% Testing m2l
ares = [0 0 0; 1 1 2; 2 3 4; 5 6 7; 8 9 10; 11 2 13];
matrix2latex(are1,'fdems_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
%% Only hosts
model = @SISV_fdem;
prev_h = @(y)(y(:,2));
num_iters = 50;
params = [betah betav gamma];
init_guess1 = params;
init_guess2 = [0.05 0.05 0.05];
[are1,pe1] = ARE(model,prev_h,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
[are2,pe2] = ARE(model,prev_h,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
save('SISV_fdems_honly_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fdems_honly_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fdems_honly_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
%% Only Vectors
model = @SISV_fdem;
prev_v = @(y)(y(:,4));
num_iters = 50;
params = [betah betav gamma];
init_guess1 = params;
init_guess2 = [0.05 0.05 0.05];
[are1,pe1] = ARE(model,prev_v,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
[are2,pe2] = ARE(model,prev_v,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
save('SISV_fdems_vonly_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fdems_vonly_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fdems_vonly_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
%% Cumulative
num_iters = 50;
cumul_b = @(y)(y(:,[5 6]));
model = @SISV_C_fdem;
y0 = [0.99 0.01 0.99 0.01 0 0];
params = [betah betav gamma];
init_guess1 = params;
init_guess2 = [0.05 0.05 0.05];
[are1,pe1] = ARE(model,cumul_b,params,y0,init_guess1,tspan,num_iters,error_levels,[],[],[],@under_hundred);
[are2,pe2] = ARE(model,cumul_b,params,y0,init_guess2,tspan,num_iters,error_levels,[],[],[],@under_hundred);
% save('SISV_fdems_c_AREs.mat','are1','are2','pe1','pe2');
matrix2latex(are1,'fdems_c_corr.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
matrix2latex(are2,'fdems_c_off.txt','rowLabels',{'0\%','1\%','5\%','10\%','20\%','30\%'},'columnLabels',{'$\beta_h$','$\beta_v$','$\gamma$'},'alignment','c','format','%.2f');
%% Scatter Plots
fig1 = Plot_Param_Estimates([pe2(4,:,1);pe2(4,:,2)],model,prev_b,[0.1,0.2],y0,tspan,@(p)([p(1) p(2) 0.1]),@OLS);
xlabel('\beta_h');
ylabel('\beta_v');
% saveas(fig1,'SISVfdem_PE_bhbv_guess_off.png');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
fig2 = Plot_Param_Estimates([pe2(4,:,1);pe2(4,:,3)],model,prev_b,[0.1,0.1],y0,tspan,@(p)([p(1) 0.2 p(2)]),@OLS);
xlabel('\beta_h');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig2,'SISVfdem_PE_bhg_guess_off.png');
fig3 = Plot_Param_Estimates([pe2(4,:,2);pe2(4,:,3)],model,prev_b,[0.2,0.1],y0,tspan,@(p)([0.1 p(1) p(2)]),@OLS);
xlabel('\beta_v');
ylabel('\gamma');
title('Fixed Populations Parameter Estimates at \sigma = 10%');
% saveas(fig3,'SISVfdem_PE_bvg_guess_off.png');
% %% Random Tests Setup
% model = @SISV_Test;
% prev_h = @prev_host;
% prev_v = @prev_vector;
% params = [0.1 0.1 1/(72.6*365) 1/46];
% y0 = [0.99 0.01 0.99 0.01];
% e_date = end_date(model,y0,params);
% tspan = linspace(0,e_date,e_date+1);
% error_levels = [0 0.01 0.05 0.1];
% [t,y] = ode45(model,tspan,y0,[],params);
% figure
% hold on
% plot(t,y(:,2));
% plot(t,y(:,4));
% hold off
% legend('I_h','I_v');
% %% Running
% init_guess = [0.15 0.15 0.00008 0.1];
% [ares,param_ests] = ARE(model,prev_h,params,y0,init_guess,tspan,200,error_levels); 
% ares
% init_guess = params;
% [ares2,~] = ARE(model,prev_h,params,y0,init_guess,tspan,200,error_levels);
% ares2
%% Error Funs
fun
function true_data = create_true_data(model,observations,params,y0,tspan)
    [~,y_true] = ode45(model,tspan,y0,[],params);
    true_data = observations(y_true);
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
    er = sum((data1 - data2).^2,'all');
end
%% Functions
function data = prev_both(y)
data = y(:,[2 4]);
end
function data = prev_host(y)
data = y(:,2);
end
function data = prev_vector(y)
data = y(:,4);
end
function date = end_date(model,y0,params)
    fulltspan = linspace(0,10000,10000);
    [t,y] = ode45(model,fulltspan,y0,[],params);
    prev_h = y(:,2);
    date = t(end);
    for i = 10:(length(t)-1)
       if ((prev_h(i+1)-prev_h(i))/(t(i+1)-t(i)) < 1e-7)
           date = t(i);
           break
       end
    end
    date = round(date);
end
function dx = SISV_fpops(t,x,z)
    dx = SISV_Model(t,x,[z(1:3) z(4) z(5) z(4) z(5)]); 
end
function dx = SISV_fdem(t,x,z)
    dx = SISV_Model(t,x,[z 0.00004 0.1 0.00004 0.1]); 
end
function dx = SISV_C_fdem(t,x,z)
    dx = SISV_Model_Cumulative(t,x,[z 0.00004 0.1 0.00004 0.1]);
end
function dx = SISV_fh(t,x,z)
    dx = SISV_Model(t,x,[z(1) z(2) z(3) 0.00004 z(4) 0.00004 z(4)]);
end
function is = under_hundred(params)
    if (params(1) > 100) || (params(2) > 100)
        is = true;
    else
        is = false;
    end
end
function dx = SISV_Test(t,x,z)
    dx = SISV_Model(t,x,[z(1) z(2) 1/10 1/(72.6*365) 1/46 z(3) z(4)]);
end