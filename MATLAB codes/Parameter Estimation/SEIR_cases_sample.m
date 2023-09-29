   %Monte Carlo Simulation from Tuncer & Le (2018) Section 4.1 
%Created: July 26, 2021
%Updated: November 30, 2021
%GLS incorporated

clear all 
close all 
clc

tic


%Number of iterations for the simulation. 
NumberofIterations=1; 

%%
%Time span
daily = 1:1:50; %Daily 
tspan=daily;



%Initial Guess
Fitted_Parameters=[4 2 1]; %beta, gamma, alpha %Peak at day 109.



%Initial Conditions   
N = 1000;
I0 = 10;
E0 = 0;
R0 = 0;
S0 = N - I0 - E0 - R0;
init_cond = [S0,E0,I0,R0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters. 

ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.  
%Fval = zeros(1,NumberofIterations); %Storage for optimization values. 

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters. 
Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels from paper.

%Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states,daily,init_cond,[],Fitted_Parameters); 


%% MC Optimization 
for IterationLevels = 1:6
IterationLevels
    
NoiseLevel = Levels(IterationLevels); %Defines noise level for current iteration. 
% use "for"  instead of "parfor" if you don't want parallel processing
for i= 1:NumberofIterations

    
Noise  = NoiseLevel*y_est(tspan(:),3)'; 
Prevalencedata = normrnd(y_est(tspan(:),3)', Noise)'; %Equation from step 2. Section 4.1. Use normal distribution. 


value=sum(Prevalencedata(:)<0); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal. 

while (value~=0)
    Prevalencedata = normrnd(y_est(tspan(:),3)', Noise)';
    value=sum(Prevalencedata(:)<0);
end

%Initial Guesses 
Initial_Guess=Fitted_Parameters; %beta, gamma, alpha. 

%Bounds for the parameters
Lowerbounds = [0 0 0];
Upperbounds=[8 4 2];

%Optimization 
%Fminsearchbnd
%options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);  
[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options); 
% value at optimal point - fval 
% whether the function - convergence rate - exitflag 

EstiParams(:,i) = EstimatedParameters'; %Stores parameters for current noise level to computer ARE.  
EstiParam(:,i, IterationLevels) = EstimatedParameters'; %Stores for each noise level
InfectData(:,i, IterationLevels) = Prevalencedata';
ExitFlag(:,i, IterationLevels) = exitflag;
Fval(:,i, IterationLevels)=fval; 

end

%Computes the ARE Score
ARE_Value = zeros(1,Number_Parameters);  %Storage for ARE calculation. 
    for i = 1:Number_Parameters
        ARE_Value(i) = (100/NumberofIterations) * sum(abs(Fitted_Parameters(i) - EstiParams(i,:)))/abs(Fitted_Parameters(i));
    end
    
    ARE(IterationLevels,:) = ARE_Value;
end

%Calcuates mean Ep, co-variances COVp, standard deviation SEp, and
%correlation coefficients rho
% Ep=mean(EstiParam');
% COVp=cov(EstiParam');
% SEp=std(EstiParam'); 
% rho=corrcoef(EstiParam');



%% Plots (Will Use Later)
% [t,y_est]=ode45(@SEIR_model_states,tspan,init_cond,[],EstimatedParameters); %Uses estimated parameters to solve system.
% 
% figure1 = figure('position', [200, 400, 1200, 400]);
% plot(tspan,data,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
% hold on 
% plot(tspan, y_est(:,3), '-b','LineWidth',3)
% title({'Fitting SEIR'},...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',15,...
% 'FontName','Times')
% ylabel({'I(t)'},...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',14,...
% 'FontName','Times')
% hold off

toc



