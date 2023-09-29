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
%% Running R0 = sqrt(2)
init_guess = [0.1 0.2 0.09996 0.00004 0.1 0.00004 0.1];
[are_c,pe_c] = ARE(model,prev_b,params,y0,init_guess,tspan,100,error_levels);
disp(are_c);
init_guess2 = [0.05 0.05 0.05 0.00004 0.1 0.00004 0.1];
[are_d,pe_d] = ARE(model,prev_b,params,y0,init_guess2,tspan,100,error_levels);
disp(are_d);
init_guess3 = [0.05 0.05 0.05 0.00002 0.05 0.00002 0.05];
[are_a,pe_a] = ARE(model,prev_b,params,y0,init_guess3,tspan,100,error_levels);
disp(are_a);
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
%% Fixing Population
model = @SISV_fpops;
params = [betah betav gamma muh muv];
init_guess1 = params;
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,100,error_levels,[],[],[]);
init_guess2 = [0.05 0.05 0.05 0.00004 0.05];
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,100,error_levels,[],[],[]);
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
params = [betah betav gamma muv];
init_guess1 = params;
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,100,error_levels,[],[],[]);
init_guess2 = [0.05 0.05 0.05 0.05];
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,100,error_levels,[],[],[]);
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
params = [betah betav gamma];
init_guess1 = params;
init_guess2 = [0.05 0.05 0.05];
[are1,pe1] = ARE(model,prev_b,params,y0,init_guess1,tspan,100,error_levels,[],[],[],[]);
[are2,pe2] = ARE(model,prev_b,params,y0,init_guess2,tspan,100,error_levels,[],[],[],[]);
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
%% R0 = 2 Setup
% We use a population of 100 hosts and 1000 vectors
% Note that Model does not change
%Start with 1/100 infected for each
model = @SISV_Model;
betah = 0.1;
betav = 0.4;
gamma = 0.09996;
muh = 0.00004;
muv = 0.1;
pih = 0.004;
piv = 100;
params = [betah betav gamma muh muv pih piv];
y0 = [99 1 990 10];
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
prev_b = @prev_both;
error_levels = [0 0.01 0.05 0.1 0.2 0.3];
%% Fixing Demographics
model = @SISV_fdem;
params = [betah betav gamma];
init_guess = [0.05 0.05 0.05];
[are,param_estimates] = ARE(model,prev_b,params,y0,init_guess,tspan,1000,error_levels,[],[],[],@isStiffd);
%% R0 = 4 Setup
% We use a population of 100 hosts and 1000 vectors
% Note that Model does not change
%Start with 1/100 infected for each
model = @SISV_Model;
betah = 0.2;
betav = 0.8;
gamma = 0.09996;
muh = 0.00004;
muv = 0.1;
pih = 0.004;
piv = 100;
params = [betah betav gamma muh muv pih piv];
y0 = [99 1 990 10];
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
prev_b = @prev_both;
error_levels = [0 0.01 0.05 0.1 0.2 0.3];
%% Fixing Demographics
model = @SISV_fdem;
params = [betah betav gamma];
init_guess = [0.05 0.05 0.05];
[are,param_estimates] = ARE(model,prev_b,params,y0,init_guess,tspan,1000,error_levels,[],[],[],@isStiffd);
%% Random Tests Setup
model = @SISV_Test;
prev_h = @prev_host;
prev_v = @prev_vector;
params = [0.1 0.1 1/(72.6*365) 1/46];
y0 = [0.99 0.01 0.99 0.01];
e_date = end_date(model,y0,params);
tspan = linspace(0,e_date,e_date+1);
error_levels = [0 0.01 0.05 0.1];
[t,y] = ode45(model,tspan,y0,[],params);
figure
hold on
plot(t,y(:,2));
plot(t,y(:,4));
hold off
legend('I_h','I_v');
%% Running
init_guess = [0.15 0.15 0.00008 0.1];
[ares,param_ests] = ARE(model,prev_h,params,y0,init_guess,tspan,200,error_levels); 
ares
init_guess = params;
[ares2,~] = ARE(model,prev_h,params,y0,init_guess,tspan,200,error_levels);
ares2
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
function is = isStiff(params)
    r0 = sqrt(params(1)*params(2)/((params(4)+params(3))*params(5)));
    is = r0 > 50;
end
function is = isStiffd(params)
is = isStiff([params 0.00004 0.1 0.00004 0.1]);
end
function dx = SISV_Test(t,x,z)
    dx = SISV_Model(t,x,[z(1) z(2) 1/10 1/(72.6*365) 1/46 z(3) z(4)]);
end