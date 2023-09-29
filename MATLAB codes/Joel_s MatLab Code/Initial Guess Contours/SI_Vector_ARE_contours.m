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
tspan = linspace(0,1000,101);
y0 = [99 1 999 1];
error_levels = 0;
num_iters = 5;
[betah_g,betav_g] = meshgrid(0.00001:0.00001:0.0002,0.0001:0.0001:0.002);
s = size(betah_g);
cols = s(2);
rows = s(1);
ares = zeros([s 7]);
%% Running
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[muh muv betah_g(r,c) betav_g(r,c) gamma pih piv],tspan,num_iters,error_levels,[],[]);
    end
end
%% Plot
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
colorbar;
%% FminCon
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[muh muv betah_g(r,c) betav_g(r,c) gamma pih piv],tspan,num_iters,error_levels,[],@fmcon);
    end
end
%% Plot
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
colorbar;
%% FminUnc
for c = 1:cols
    for r = 1:rows
        ares(r,c,:) = ARE(model,prev,params,y0,[muh muv betah_g(r,c) betav_g(r,c) gamma pih piv],tspan,num_iters,error_levels,[],@fmunc);
    end
end
%% Plot
figure
contourf(betah_g,betav_g,ares(:,:,3),[0 0.1 0.5 1 2 5 10])
colorbar;
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
function newpars = fmcon(error_fun,init_guess)
    options = optimoptions('fmincon','Display','none');
    newpars = fmincon(error_fun,init_guess,[],[],[],[],[],[],[],options);%For some reason, giving it 0 as a lower bound makes it significantly worse
end
function newpars = fmunc(error_fun,init_guess)
    options = optimoptions('fminunc','Display','none');
    newpars = fminunc(error_fun,init_guess,options);
end