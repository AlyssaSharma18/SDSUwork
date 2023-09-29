%% Model settings
model = @SIVector_Model;
prev = @prev_both;
pih = 0.004;
piv = 90;
muh = .00004;
muv = 0.09;
betah = 0.0001;
betav = 0.001;
gamma = .08996;
params = [muh muv betah betav gamma pih piv];
init_guess = params;
tspan = linspace(0,1000,101);
y0 = [99 1 999 1];
error_levels = [0 0.01 0.05 0.1 0.2 0.3];
%% Both prevalence
num_iters = 10;
are = ARE(model,prev,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
%writematrix(are,'SIV_AREs_default.xlsx');
%% Only Host Prevalence
num_iters = 10;
are = ARE(model,@prev_host,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
%writematrix(are,'SIV_AREs_honly_def.xlsx');
%% Only Vector Prevalence
num_iters = 10;
are = ARE(model,@prev_vec,params,y0,init_guess,tspan,num_iters,error_levels,[],[]);
%writematrix(are,'SIV_AREs_vonly_def.xlsx');
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