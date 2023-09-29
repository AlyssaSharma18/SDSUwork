%Practice comparing SIR vs SEIR 

t=1:1:200;
paras = [1/3000,1/6,1/2,0];
parase = [1/3000,1/6,1/2,1];
init = [999,0,1,0];

[t1,sirs]=ode45(@compartmentsODEs,t,init,[],paras);
[t2,seirs] = ode45(@compartmentsODEs,t, init, [],parase);
%%
T=100;
figure
plot(t(1:T),sirs(1:T,3),'b','LineWidth',1.5)
hold on
plot(t(1:T),seirs(1:T,3),'r','LineWidth',1.5)

xlabel('Time')
ylabel('Infected')
legend('SIR','SEIR')
title([' parameters: \beta = ',num2str(paras(1)),', \alpha = ',num2str(paras(2)),', k = ', num2str(paras(3))])
hold off