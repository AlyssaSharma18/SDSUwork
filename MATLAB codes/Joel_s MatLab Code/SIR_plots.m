%list_params = [0.01,1;0.1,1;0.25,1;0.5,1;0.95,1;0.99,1];
%Have it also plot the top?
list_params = [0.5,1];
figure
hold on
for params = list_params'
    num_points = 1000;
    tspan = linspace(1,100,num_points);
    y0 = [.999 .001 0];
    [t,y] = ode45(@(t,x) SIR_model(t,x,params),tspan,y0);
    alpha = params(1);
%     if (alpha > 0.99)
%         prev = 1000 .* y(:,2);
%     else
    prev = y(:,2);
%     end
    plot(tspan,prev)
end
% legend('0.01','0.1','0.5','1','2','3.48(*1000)');
% title('Infected by \alpha value at \beta = 3.5');
% ylabel('Infections');
% xlabel('Day');
hold off
%saveas(gcf,'peak_infection_time_plots.png')