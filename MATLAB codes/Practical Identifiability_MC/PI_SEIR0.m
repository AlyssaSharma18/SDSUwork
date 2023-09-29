clear all
%Number of iterations for the simulation.
NumberofIterations=10;

%%
%Time span
T=50;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.0001 0.2 0.03]; %beta, gamma, alpha %Peak at day 109.

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

    NoiseLevel = Levels(IterationLevels); %Defines noise level for current iteration.
    % use "for"  instead of "parfor" if you don't want parallel processing
    for i= 1:NumberofIterations
    outp=sprintf('it is current at noise level %d, iteration number %d',IterationLevels,i);
    disp(outp)
    % these two lines are used for error check when you run a small number
    % of iterations under each noise level. When you run 1000 simulations,
    % better not have them, or only print out every 100 iterations 

        Noise  = NoiseLevel*y_est(tspan(:),3)';
        Prevalencedata = normrnd(y_est(tspan(:),3)', Noise)'; %Equation from step 2. Section 4.1. Use normal distribution.


        value=sum(Prevalencedata(:)<-0.001); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal.

        while (value~=0)
            Prevalencedata = normrnd(y_est(tspan(:),3)', Noise)';
            value=sum(Prevalencedata(:)<0);
        end

        %Initial Guesses
        Initial_Guess=Fitted_Parameters; %beta, gamma, alpha.

        %Bounds for the parameters
        Lowerbounds = [0 0 0];
        Upperbounds=[1 1 1];

        %Optimization
        %Fminsearchbnd
        %[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options);
        options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
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

%%  SEIR model
function [dx]  = SEIR_model_states(t,x,z)

%parameters to be estimated
beta = z(1);
gamma = z(2);
alpha = z(3);


S = x(1); E = x(2); I = x(3); R= x(4); N=1000;

%IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0


dx = zeros(4,1);
%dS
dx(1) = -beta*S*I;
%dE
dx(2) = beta*S*I-gamma*E;
%dI
dx(3) = gamma*E-alpha*I;
%dR
dx(4) = alpha*I;
end


%%

function  error_in_data = err_in_dataSEIR_GLS(par,data,T)
%This is the error function for prevalence.
%Initial Conditions
N = 1000;
I0 = 10;
E0 = 0;
R0 = 0;
S0 = N - I0 - E0 - R0;
init_cond = [S0,E0,I0,R0];



daily = 1:1:T; %Daily


tspan=daily;

[t,y]=ode45(@SEIR_model_states,daily,init_cond,[],par);

Infected=y(tspan(:),3);
weight=1./Infected.^2;
error_in_data = sum(weight.*(Infected-data).^2);

end






