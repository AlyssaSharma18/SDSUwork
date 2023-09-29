clear all
%Number of iterations for the simulation.
NumberofIterations=2;

%%
%Time span
T=365;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.00019 0.3162 0.3162 0.05 0.05]; %betaaa, betacc, betaac, betaca, gammac, gammaa, epsilonc, epsilona
%Peak at day 109.

%Initial Conditions
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
init_cond = [Sc0, Sa0, Ec0, Ea0, Ic0, Ia0, Rc0, Ra0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters.

ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.
%Fval = zeros(1,NumberofIterations); %Storage for optimization values.

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters.
Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels from paper.

%Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states_AdultChild_a,daily,init_cond,[],Fitted_Parameters);


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
        for j= (NumberofIterations+1): (NumberofIterations*2)
        Noisec  = NoiseLevel*y_est(tspan(:),5)';
        Noisea = NoiseLevel*y_est(tspan(:),6)';
        Prevalencedatac = normrnd(y_est(tspan(:),5)', Noisec)'; %Equation from step 2. Section 4.1. Use normal distribution.
        Prevalencedataa = normrnd(y_est(tspan(:),6)', Noisea)';


        valuec=sum(Prevalencedatac(:)<-0.001); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal.
        valuea=sum(Prevalencedataa(:)<-0.001);

        while (valuec~=0)
            Prevalencedatac = normrnd(y_est(tspan(:),5)', Noisec)';
            valuec=sum(Prevalencedatac(:)<-0.001);
        end

        while (valuea~=0)
            Prevalencedataa = normrnd(y_est(tspan(:),6)', Noisea)';
            valuea=sum(Prevalencedataa(:)<-0.001);
        end

        %Initial Guesses
        Initial_Guess=Fitted_Parameters; %beta, gamma, alpha.

        %Bounds for the parameters
        Lowerbounds = [0 0 0 0 0];
        Upperbounds=[1 1 1 1 1];

        %Optimization
        %Fminsearchbnd
        %[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options);
        options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
        [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_dataSEIR_GLS_AdultChild(Initial_Guess,Prevalencedatac,Prevalencedataa, T),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);


        EstiParams(:,i) = EstimatedParameters'; %Stores parameters for current noise level to computer ARE.
        EstiParam(:,i, IterationLevels) = EstimatedParameters'; %Stores for each noise level
        InfectData(:,i, IterationLevels) = Prevalencedatac';
        InfectData(:,j, IterationLevels) = Prevalencedataa';
        ExitFlag(:,i, IterationLevels) = exitflag;
        Fval(:,i, IterationLevels)=fval;
        end
    end

    %Computes the ARE Score
    ARE_Value = zeros(1,Number_Parameters);  %Storage for ARE calculation.
    for i = 1:Number_Parameters
        ARE_Value(i) = (100/NumberofIterations) * sum(abs(Fitted_Parameters(i) - EstiParams(i,:)))/abs(Fitted_Parameters(i));
    end

    ARE(IterationLevels,:) = ARE_Value;
end

%%  SEIR model
function [dx]  = SEIR_model_states_AdultChild_a(t,x,z)

    %parameters to be estimated
    beta = z(1); %0.00019
    betaaa = beta;
    betacc = (3/2)*250*beta;
    betaac = (1/3)*750*beta;
    betaca = (3/4)*beta;
    gammaa = z(2); %0.3162; %z(5);
    gammac = z(3); %0.3162; %z(6);
    epsilonc = z(4); %0.05; %z(7);
    epsilona = z(5); %0.05; %z(8);
    mu = 0.00003;
    f = mu/6000;
    Nc = 250;
    Na = 750;
      
    Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8);
    
   %IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0
   
   
    dx = zeros(8,1);
    %dSc
    dx(1) = mu*Na - (betacc*Ic/Nc + betaac*Ia/Na)*Sc - f*Sc;
    %dSa 
    dx(2) = f*Sc - (betaaa*Ia+betaca*Ic)*Sa - mu*Sa;
    %dEc
    dx(3) = (betacc*Ic/Nc + betaac*Ia/Na)*Sc - epsilonc*Ec - f*Ec;
    %dEa
    dx(4) = (betaaa*Ia + betaca*Ic)*Sa - epsilona*Ea + f*Ec - mu*Ea;
    %dIc
    dx(5) = epsilonc*Ec - f*Ic - gammac*Ic;
    %dIa 
    dx(6) = f*Ic+ epsilona*Ea - gammaa*Ia - mu*Ia;
    %dRc
    dx(7) = gammac*Ic - f*Rc;
    %dRa
    dx(8) = gammaa*Ia +f*Rc - mu*Ra;

    
end


%%

function  error_in_data = err_in_dataSEIR_GLS_AdultChild(par, datac, dataa, T)
%This is the error function for prevalence.    
%Initial Conditions
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
init_cond = [Sc0, Sa0, Ec0, Ea0, Ic0, Ia0, Rc0, Ra0];



daily = 1:1:365; %Daily  


tspan=daily; 

[t,y]=ode45(@SEIR_model_states_AdultChild_a, daily, init_cond,[],par); 

Infectedc=y(tspan(:),5);
Infecteda=y(tspan(:),6);
weightc=1./Infectedc.^2;  
weighta=1./Infecteda.^2; 
error_in_data = sum(weightc.*(Infectedc-datac).^2) + sum(weighta.*(Infecteda-dataa).^2);
end