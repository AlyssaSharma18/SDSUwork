clear all
%Number of iterations for the simulation.
NumberofIterations=100;

%%
%Time span
T=50;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.3162 0.3162 0.5 0.5]; % beta_aa beta_cc beta_ac beta_ca gamma_c gamma_a epsilon_c epsilon_a    %Peak at day 109.

%Initial Conditions
Sc0 = 240;
Sa0 = 740;
Ec0 = 0;
Ea0 = 0;
Ic0 = 10;
Ia0 = 10;
Rc0 = 0;
Ra0 = 0;
Nc = Sc0 + Ic0 + Ec0 + Rc0;
Na = Sa0 + Ia0 + Ea0 + Ra0;
init_cond = [Sc0,Sa0,Ec0,Ea0,Ic0,Ia0,Rc0,Ra0];
% %N = 1000;
% I0 = 10;
% E0 = 0;
% R0 = 0;
% S0 = N - I0 - E0 - R0;
% init_cond = [S0,E0,I0,R0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters.

ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.
%Fval = zeros(1,NumberofIterations); %Storage for optimization values.

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters.
Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels from paper.

%Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states_AdultChild,daily,init_cond,[],Fitted_Parameters);


%% MC Optimization
for IterationLevels = 1:6

    NoiseLevel = Levels(IterationLevels); %Defines noise level for current iteration.
    % use "for"  instead of "parfor" if you don't want parallel processing
    for i= 1:NumberofIterations
    outp=sprintf('it is current at noise level %d, iteration number %d',IterationLevels,i);
    disp(outp)
    % these two lines are used for error check when you run a small number
    % of iterations under each noise level. When you run 1000 simulations,
    % better not have them, or only print out every 100 iterations 

        Noise  = NoiseLevel*y_est(tspan(:),5)';
        Prevalencedata = normrnd(y_est(tspan(:),5)', Noise)'; %Equation from step 2. Section 4.1. Use normal distribution.


        value=sum(Prevalencedata(:)<-0.001); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal.

        while (value~=0)
            Prevalencedata = normrnd(y_est(tspan(:),5)', Noise)';
            value=sum(Prevalencedata(:)<0);
        end

        %Initial Guesses
        Initial_Guess=Fitted_Parameters; %beta, gamma, alpha.

        %Bounds for the parameters
        Lowerbounds = [0 0 0 0];
        Upperbounds=[1 1 1 1];

        %Optimization
        %Fminsearchbnd
         options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
        %[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options);
       
        [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata,T),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);


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


%% SEIR Adult Child
function [dx] = SEIR_model_states_AdultChild(t,x,z)

%parameters to be estimated

betaaa = 0.0019;%z(1);
betacc = 0.18; %z(2);
betaac = 0.013; %z(3);
betaca = 0.00013; %z(4);
gammac = z(1);
gammaa = z(2);
epsilonc = z(3);
epsilona = z(4);
mu = 0.00003;
f = mu/6000;
Nc = 250;
Na = 750;


Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8);

dx = zeros(8,1);
%dSc
dx(1) = mu*Na - (betacc*Ic/Nc + betaac*Ia/Na)*Sc - f*Sc;
%dSa
dx(2) = f*Sc - (betaaa*Ia + betaca*Ic)*Sa - mu*Sa;
%dEc
dx(3) = (betacc*Ic/Nc + betaac*Ia/Na)*Sc - epsilonc*Ec - f*Ec;
%dEa
dx(4) = (betaaa*Ia + betaca*Ic)*Sa - epsilona*Ea + f*Ec - mu*Ea;
%dIc 
dx(5) = epsilonc*Ec - f*Ic - gammac*Ic;
%dIa 
dx(6) = f*Ic + epsilona*Ea - gammaa*Ia - mu*Ia;
%dRc
dx(7) = gammac*Ic - f*Rc;
%dRa
dx(8) = gammaa*Ia + f*Rc - mu*Ra;

end


%%

function  error_in_data = err_in_dataSEIR_GLS(par,data,T)
%This is the error function for prevalence.
%Initial Conditions
% N = 1000;
% I0 = 10;
% E0 = 0;
% R0 = 0;
% S0 = N - I0 - E0 - R0;
% init_cond = [S0,E0,I0,R0];

Sc0 = 240;
Sa0 = 740;
Ic0 = 10;
Ia0 = 10;
Ec0 = 0;
Ea0 = 0;
Rc0 = 0;
Ra0 = 0;
Nc = Sc0 + Ic0 + Ec0 + Rc0;
Na = Sa0 + Ia0 + Ea0 + Ra0;
init_cond = [Sc0,Sa0,Ec0,Ea0,Ic0,Ia0,Rc0,Ra0];



daily = 1:1:T; %Daily


tspan=daily;

[t,y]=ode45(@SEIR_model_states_AdultChild,daily,init_cond,[],par);

Infected=y(tspan(:),5); %6 for adults
weight=1./Infected.^2;
error_in_data = sum(weight.*(Infected-data).^2);

end





