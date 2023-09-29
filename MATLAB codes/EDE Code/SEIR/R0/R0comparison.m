%% Main Value
clear all;

beta = 0.00027; 
betaaa = beta; 
betacc = (1/4)*(24515/42739)*250*beta;
betaac = (1/6)*750*beta;
betaca = (1/3)*(24515/42739)*beta;
% Use with seirR0.mlx code
trueR0 = seirR0([betaaa betacc betaac betaca 0.074 0.1 0.3 0.2]); % This calculates the R0 for our model at the true parameters, given a dependent beta

%% 

changebeta = linspace(-0.004,0.004,75);
%changebeta = [-0.00005 -0.00004 -0.00003 -0.00002 -0.00001 0.0000 0.00001 0.00002 0.00003 0.00004]; % This is how you vary your beta values, beta plus or minus something
for i = 1:75
    betaaa(i) = beta + changebeta(i); %we increase betaaa and decrease betaca by the same step
    betaca(i) = (1/3)*(24515/42739)*(beta - changebeta(i));
    pointwiseR0(i) = seirR0([betaaa(i) betacc betaac betaca(i) 0.074 0.1 0.3 0.2]); %calculate the new R0 value for the changed parameters
    if pointwiseR0(i) < trueR0 %our color function returns 3 different values depending on if the new R0 is greater than, equal to, or less than the true R0
        c(i) = -1;
    elseif pointwiseR0(i) == trueR0
        c(i) = 0;
    else
        c(i) = 1;
    end

end

scatter(betaaa,betaca,50,c,'filled');

% %% Calculate R0
% 
% function R0 = seirR0(par)
% betaaa = par(1);
% betacc = par(2);
% betaac = par(3);
% betaca = par(4);
% gammaa = par(5);
% gammac = par(6); 
% epsilonc = par(7);
% epsilona = par(8);
% mu = 0.00007;
% f = 0.0005;
% Nc = 250;
% Na = 750;
% 
% eig1 = betacc*(epsilonc/((epsilonc+f)*(f+gammac))) + Nc*betaac*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/(Na*(epsilonc+f)*(f+gammac)*(epsilona+mu)*(gammaa+mu)));
% eig2 = Nc*betaac*(epsilona/(Na*(epsilona + mu)*(gammaa + mu)));
% eig3 = Na*betaca.*(epsilonc/((epsilonc+f)*(f+gammac))) + Na*betaaa*((epsilona*f^2 + epsilona*epsilonc*f + epsilona*f*gammac + epsilonc*f*mu)/((epsilonc + f)*(f + gammac)*(epsilona + mu)*(gammaa + mu)));
% eig4 = Na*betaaa*(epsilona/((epsilona+mu)*(gammaa+mu)));
% 
% R0 = eig1/2 + eig4/2 + ((eig1*eig1 + eig4*eig4 - 2*eig1*eig4 + 4*eig2*eig3)^(1/2));
% end
