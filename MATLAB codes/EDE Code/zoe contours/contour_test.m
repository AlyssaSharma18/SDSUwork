Nc = 250;
Na = 750;
beta = 0.00027;
betaca = beta; 
betaaa = beta; %beta - delta2
betacc = beta*Nc;
betaac = beta*Na;
epsilonc = 0.25;
epsilona = 0.25; %average
gammac = 0.087;
gammaa = 0.087;
mu = 0.00007;
f = 0.0005;

delta1 = linspace(0, 2*beta,51);
delta2 = linspace(0, 2*beta,51);
[A,B] = meshgrid(delta1,delta2); %betaca, betaaa

eig1 = betacc*(epsilonc/((epsilonc+f)*(f+gammac))) + Nc*betaac*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/(Na*(epsilonc+f)*(f+gammac)*(epsilona+mu)*(gammaa+mu)));
eig2 = Nc*betaac*(epsilona/(Na*(epsilona + mu)*(gammaa + mu)));
eig3 = Na*A*(epsilonc/((epsilonc+f)*(f+gammac))) + Na*B*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/((epsilonc + f)*(f + gammac)*(epsilona + mu)*(gammaa + mu)));
eig4 = Na*B*(epsilona/((epsilona+mu)*(gammaa+mu)));
 
Ro = eig1/2 + eig4/2 + ((eig1.*eig1 + eig4.*eig4 - 2*eig1.*eig4 + 4*eig2.*eig3).^(1/2))/2;

figure()
hold on
contourf(A,B,Ro)
title("\beta_{ca} vs \beta_{aa} with N_c = 250 and N_a = 750")
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")
plot(beta,beta,'r*')
colorbar
hold off
%%
for i = 1:numel(Ro)
    if Ro(i) > Ro(26,26)
        Ro(i) = 1;

    elseif Ro(i) <= Ro(26,26) && Ro(i) >= 1
        Ro(i) = 0;

    elseif Ro(i) < 1
        Ro(i) = -1;
    end 
end 

%%
figure()
hold on
%heatmap(delta1,delta2,Ro)
contourf(A,B,Ro)
plot(beta,beta,'r*')
title(["\beta_{ca} vs \beta_{aa} with N_c = 250 and N_a = 750",  "In comparision with simpler SEIR model"])
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")
colorbar
hold off

%%
figure()
hM = heatmap(delta1, delta2, Ro, 'Colormap',jet)
hM.NodeChildren(3).YDir='normal';                   % turn Y-Axis normal direction
title(["\beta_{ca} vs \beta_{aa} with N_c = 250 and N_a = 750",  "In comparision with simpler SEIR model"])
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")

%%
Nc = 500;
Na = 500;

delta1 = linspace(0, 2*beta,51);
delta2 = linspace(0, 2*beta,51);
[A,B] = meshgrid(delta1,delta2); %betaca, betaaa

eig1 = betacc*(epsilonc/((epsilonc+f)*(f+gammac))) + Nc*betaac*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/(Na*(epsilonc+f)*(f+gammac)*(epsilona+mu)*(gammaa+mu)));
eig2 = Nc*betaac*(epsilona/(Na*(epsilona + mu)*(gammaa + mu)));
eig3 = Na*A*(epsilonc/((epsilonc+f)*(f+gammac))) + Na*B*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/((epsilonc + f)*(f + gammac)*(epsilona + mu)*(gammaa + mu)));
eig4 = Na*B*(epsilona/((epsilona+mu)*(gammaa+mu)));
 
Ro = eig1/2 + eig4/2 + ((eig1.*eig1 + eig4.*eig4 - 2*eig1.*eig4 + 4*eig2.*eig3).^(1/2))/2;

figure()
hold on
contourf(A,B,Ro)
title("\beta_{ca} vs \beta_{aa} with N_c = N_a = 500")
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")
plot(beta,beta,'r*')
colorbar
hold off


%%
for i = 1:numel(Ro)
    if Ro(i) > Ro(26,26)
        Ro(i) = 1;

    elseif Ro(i) <= Ro(26,26) && Ro(i) >= 1
        Ro(i) = 0;

    elseif Ro(i) < 1
        Ro(i) = -1;
    end 
end 

%%
figure()
hold on
%heatmap(delta1,delta2,Ro)
contourf(A,B,Ro)
plot(beta,beta,'r*')
title(["\beta_{ca} vs \beta_{aa} with N_c = N_a = 500",  "In comparision with simpler SEIR model"])
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")
colorbar
hold off

%%
figure()
hM = heatmap(delta1, delta2, Ro, 'Colormap',jet)
hM.NodeChildren(3).YDir='normal';                   % turn Y-Axis normal direction
title(["\beta_{ca} vs \beta_{aa} with N_c = N_a = 500",  "In comparision with simpler SEIR model"])
xlabel("\beta_{ca}")
ylabel("\beta_{aa}")




