figure
hold on
scatter(EstiParam(4,:,2),EstiParam(5,:,2))
scatter(Fitted_Parameters(4), Fitted_Parameters(5), 'filled')
xline(Fitted_Parameters(4)- 0.01*Fitted_Parameters(4))
xline(Fitted_Parameters(4)+ 0.01*Fitted_Parameters(4))
yline(Fitted_Parameters(5)- 0.01*Fitted_Parameters(5))
yline(Fitted_Parameters(5)+ 0.01*Fitted_Parameters(5))
title('\epsilon_{C} and \epsilon_{A} \sigma=1')
xlabel('\epsilon_{C}')
ylabel('\epsilon_{A}')
axis([0.29 .31 0.195 .205])
hold off

figure
hold on
scatter(EstiParam(4,:,3),EstiParam(5,:,3))
scatter(Fitted_Parameters(4), Fitted_Parameters(5), 'filled')
xline(Fitted_Parameters(4)- 0.05*Fitted_Parameters(4))
xline(Fitted_Parameters(4)+ 0.05*Fitted_Parameters(4))
yline(Fitted_Parameters(5)- 0.05*Fitted_Parameters(5))
yline(Fitted_Parameters(5)+ 0.05*Fitted_Parameters(5))
title('\epsilon_{C} and \epsilon_{A} \sigma=5')
xlabel('\epsilon_{C}')
ylabel('\epsilon_{A}')
axis([0.25 .35 0.17 .23])
hold off

figure
hold on
scatter(EstiParam(4,:,4),EstiParam(5,:,4))
scatter(Fitted_Parameters(4), Fitted_Parameters(5), 'filled')
xline(Fitted_Parameters(4)- 0.1*Fitted_Parameters(4))
xline(Fitted_Parameters(4)+ 0.1*Fitted_Parameters(4))
yline(Fitted_Parameters(5)- 0.1*Fitted_Parameters(5))
yline(Fitted_Parameters(5)+ 0.1*Fitted_Parameters(5))
title('\epsilon_{C} and \epsilon_{A} \sigma=10')
xlabel('\epsilon_{C}')
ylabel('\epsilon_{A}')
axis([0.2 .4 0.15 .25])
hold off

figure
hold on
scatter(EstiParam(4,:,5),EstiParam(5,:,5))
scatter(Fitted_Parameters(4), Fitted_Parameters(5), 'filled')
xline(Fitted_Parameters(4)- 0.2*Fitted_Parameters(4))
xline(Fitted_Parameters(4)+ 0.2*Fitted_Parameters(4))
yline(Fitted_Parameters(5)- 0.2*Fitted_Parameters(5))
yline(Fitted_Parameters(5)+ 0.2*Fitted_Parameters(5))
title('\epsilon_{C} and \epsilon_{A} \sigma=20')
xlabel('\epsilon_{C}')
ylabel('\epsilon_{A}')
axis([0.15 .6 0.1 .3])
hold off

figure
hold on
scatter(EstiParam(4,:,6),EstiParam(5,:,6))
scatter(Fitted_Parameters(4), Fitted_Parameters(5), 'filled')
xline(Fitted_Parameters(4)- 0.3*Fitted_Parameters(4))
xline(Fitted_Parameters(4)+ 0.3*Fitted_Parameters(4))
yline(Fitted_Parameters(5)- 0.3*Fitted_Parameters(5))
yline(Fitted_Parameters(5)+ 0.3*Fitted_Parameters(5))
title('\epsilon_{C} and \epsilon_{A} \sigma=30')
xlabel('\epsilon_{C}')
ylabel('\epsilon_{A}')
axis([0.1 1.2 0.1 .3])
hold off