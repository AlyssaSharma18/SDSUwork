%% R0
[alpha,beta] = meshgrid(0.05:0.01:1,0.01:0.01:8);
R0 = beta./alpha;
contourf(alpha,beta,R0)
%Try for seir
%% Redoing
[alpha,beta] = meshgrid(0.05:0.1:4,0.05:0.1:8);
s = size(alpha);
y = zeros([s 1000 3]);
cols = s(2);
rows = s(1);
parfor c = 1:cols
    for r = 1:rows
        y(r,c,:,:) = data(alpha(1,c),beta(r,1),0.5);
    end
end
I_max = zeros(s);
t_max = zeros(s);
s_inf = zeros(s);
for c = 1:cols
    for r = 1:rows
        [I_max(r,c),t_max(r,c)] = max(y(r,c,:,2));
        s_inf(r,c) = y(r,c,end,1);
    end
end
%% Graph I_max
figure
contourf(alpha,beta,I_max,[0 0.011 0.05 0.1 0.2 0.3 0.5 0.6 0.7 0.8 0.9])
title('Peak Infections');
xlabel('\alpha');
ylabel('\beta');
colorbar;
%saveas(gcf,'SIR_Peak_Infections.png')
%% Graph T_max
figure
contourf(alpha,beta,t_max,[0 1.1,2.1,3.1 4.1 5.1,10.1,20.1,30.1,50.1])
title('Time of Peak Infection')
xlabel('\alpha')
ylabel('\beta')
colorbar;
%saveas(gcf,'SIR_Peak_Infect_Time.png')
%% Graph S_inf
figure
contourf(alpha,beta,s_inf,[0 0.01 0.02 0.1 0.2 0.35 0.5 0.7 0.8 0.9 0.95 0.98 1])
title('Number Susceptible After 365 Days')
xlabel('\alpha')
ylabel('\beta')
colorbar;
%saveas(gcf,'SIR_S_inf.png')
function y = data(alpha,beta)
   tspan = linspace(1,100,1000);
    y0 = [.999 .001 0];
   [~,y] = ode45(@SIR_model,tspan,y0,[],[alpha,beta]);
end
function y = SEIR_data(alpha,beta,gamma)
    tspan = linspace(1,100,1000);
    y0 = [.999, 0, .001, 0];
    [~,y] = ode45(@SEIR_model,tspan,y0,[],[alpha,beta,gamma]);
end
