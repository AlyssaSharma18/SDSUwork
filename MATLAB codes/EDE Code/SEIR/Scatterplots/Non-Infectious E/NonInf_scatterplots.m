figure
hold on
scatter(NonInf_EstiParam(1,:,2),NonInf_EstiParam(2,:,2))
scatter(NonInf_FitParam(1), NonInf_FitParam(2), 'filled')
xline(NonInf_FitParam(1)- 0.01*NonInf_FitParam(1))
xline(NonInf_FitParam(1)+ 0.01*NonInf_FitParam(1))
yline(NonInf_FitParam(2)- 0.01*NonInf_FitParam(2))
yline(NonInf_FitParam(2)+ 0.01*NonInf_FitParam(2))
title('\beta_{AA} and \beta_{CC} \sigma=1')
xlabel('\beta_{AA}')
ylabel('\beta_{CC}')
axis([0.000267 .000273 0.00985 .01015])
hold off

figure
hold on
scatter(NonInf_EstiParam(1,:,3),NonInf_EstiParam(2,:,3))
scatter(NonInf_FitParam(1), NonInf_FitParam(2), 'filled')
xline(NonInf_FitParam(1)- 0.05*NonInf_FitParam(1))
xline(NonInf_FitParam(1)+ 0.05*NonInf_FitParam(1))
yline(NonInf_FitParam(2)- 0.05*NonInf_FitParam(2))
yline(NonInf_FitParam(2)+ 0.05*NonInf_FitParam(2))
title('\beta_{AA} and \beta_{CC} \sigma=5')
xlabel('\beta_{AA}')
ylabel('\beta_{CC}')
%axis([0.25 .35 0.17 .23])
hold off

figure
hold on
scatter(NonInf_EstiParam(1,:,4),NonInf_EstiParam(2,:,4))
scatter(NonInf_FitParam(1), NonInf_FitParam(2), 'filled')
xline(NonInf_FitParam(1)- 0.1*NonInf_FitParam(1))
xline(NonInf_FitParam(1)+ 0.1*NonInf_FitParam(1))
yline(NonInf_FitParam(2)- 0.1*NonInf_FitParam(2))
yline(NonInf_FitParam(2)+ 0.1*NonInf_FitParam(2))
title('\beta_{AA} and \beta_{CC} \sigma=10')
xlabel('\beta_{AA}')
ylabel('\beta_{CC}')
%axis([0.2 .4 0.15 .25])
hold off

figure
hold on
scatter(NonInf_EstiParam(1,:,5),NonInf_EstiParam(2,:,5))
scatter(NonInf_FitParam(1), NonInf_FitParam(2), 'filled')
xline(NonInf_FitParam(1)- 0.2*NonInf_FitParam(1))
xline(NonInf_FitParam(1)+ 0.2*NonInf_FitParam(1))
yline(NonInf_FitParam(2)- 0.2*NonInf_FitParam(2))
yline(NonInf_FitParam(2)+ 0.2*NonInf_FitParam(2))
title('\beta_{AA} and \beta_{CC} \sigma=20')
xlabel('\beta_{AA}')
ylabel('\beta_{CC}')
%axis([0.15 .6 0.1 .3])
hold off

figure
hold on
scatter(NonInf_EstiParam(1,:,6),NonInf_EstiParam(2,:,6))
scatter(NonInf_FitParam(1), NonInf_FitParam(2), 'filled')
xline(NonInf_FitParam(1)- 0.3*NonInf_FitParam(1))
xline(NonInf_FitParam(1)+ 0.3*NonInf_FitParam(1))
yline(NonInf_FitParam(2)- 0.3*NonInf_FitParam(2))
yline(NonInf_FitParam(2)+ 0.3*NonInf_FitParam(2))
title('\beta_{AA} and \beta_{CC} \sigma=30')
xlabel('\beta_{AA}')
ylabel('\beta_{CC}')
%axis([0.1 1.2 0.1 .3])
hold off