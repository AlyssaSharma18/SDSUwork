%%Plots host and vector infections for a list of parameter lists
%list_params = [0.0001 0.1 1 0.01 0.1;0.0001 0.1 1 0.05 0.1;0.0001 0.1 0.1 0.01 0.1];
pih = 0.01;%Host birth rate
piv = 100;%Vector birth rate
muh = 0.0001;%Host death rate
muv = 0.1;%Vector death rate
betah = .00025;%Host infection rate
betav = 0.01;%Vector infection rate
gamma = 0.1;%Host recovery rate
list_params = [muh muv betah betav gamma pih piv];
%R0 calculation
r0 = sqrt(betah*betav*pih*piv/((muh+gamma)*(muv)*muv*muh))
Nh = pih/muh
Nv = piv/muv
iheq =betah*Nh*Nv*(1-1/(r0^2))/(muh+gamma+betah*Nv)
iveq = betav*Nh*Nv*(1-1/(r0^2))/(muv+betav*Nh)
figure
hold on
for params = list_params'
    %Solving
    num_points = 10000;
    tspan = linspace(1,100,num_points);
    y0 = [100 0 75.9 .1];
    [t,y] = ode45(@(t,x) SIVector_Model(t,x,params),tspan,y0);
    %Plotting
    plot(tspan,y(:,[2 4]),'LineWidth',1.5)
end
%Adding style to graph
%legend('R0 = 1:H','R0 = 1:V','R0 = 2.22:H','R0 = 2.22:V','R0 = .31:H','R0 = .31:V');
title('Behaviours for the SIS Vector Model');
ylabel('Infections');
xlabel('Day');
% axis([1 200 0 1]);
legend('I_h','I_v')
hold off
% saveas(gcf,'vector_SIS_plots.png')