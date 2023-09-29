clear all
%Number of iterations for the simulation.
NumberofIterations=10;
%%
%Time span
T=100;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.5 0.5 0.5 0.5]; %[0.1 0.074 0.3 0.2];
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
Nc0 = Sc0 + Ic0 + Ec0 + Rc0;
Na0 = Sa0 + Ia0 + Ea0 + Ra0;
init_cond = [Sc0, Sa0, Ec0, Ea0, Ic0, Ia0, Rc0, Ra0, Nc0, Na0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters.

ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.
%Fval = zeros(1,NumberofIterations); %Storage for optimization values.

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters.

%Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels from paper.
Levels = [0.3];

%Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states_AdultChild_infectE,daily,init_cond,[],Fitted_Parameters);

%% MC Optimization
for IterationLevels = 1
    NoiseLevel = Levels(IterationLevels); %Defines noise level for current iteration.
    % use "for"  instead of "parfor" if you don't want parallel processing
    for i= 1:NumberofIterations
    outp=sprintf('it is current at noise level %d, iteration number %d',IterationLevels,i);
    disp(outp)
    % these two lines are used for error check when you run a small number
    % of iterations under each noise level. When you run 1000 simulations,
    % better not have them, or only print out every 100 iterations 
        for j= (NumberofIterations+1): (NumberofIterations*2)
       
        Prevalencedatac = y_est(:,5) + NoiseLevel*(y_est(:,5));
        Prevalencedataa = y_est(:,6) + NoiseLevel*(y_est(:,6));
        
        %Initial Guesses
        Initial_Guess=Fitted_Parameters; %beta, gamma, alpha.

        %Bounds for the parameters
        Lowerbounds = [0 0 0 0];
        Upperbounds=[2 2 2 2];

        %Optimization
        %Fminsearchbnd
        %[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options);
        options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
        [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_dataSEIR_GLS_AdultChild_infectE(Initial_Guess,Prevalencedatac,Prevalencedataa, T, NoiseLevel),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);


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
function [dx]  = SEIR_model_states_AdultChild_infectE(t,x,z)

    %parameters to be estimated
    betaaa =0.00027; %5*10^-7 to 2*10^-5
    betacc = 0.01; %0.1 to 1
    betaac = 0.01; %.05 to 1 <betacc
    betaca = 0.00005; % 1*10^-8 to 5*10^-7
    gammac = z(1); %0.1; % 0.1 %.02 to 1
    gammaa = z(2); %0.074; % 0.074 % .01 to 1
    epsilonc = z(3); %0.3;% 0.3 % .04 to 1
    epsilona = z(4); %0.2;% 0.2 % .02 to 1
    xicc = 0.0067;
    xiaa = 0.00018;
    xica = 0.000033;
    xiac = 0.0067;
    mu = 0.00008; %4*10^-5 to 8*10^-5
    f = 0.0005; %1*10^-4 to 5*10^-4
    Nc0 = 250;
    Na0 = 750;

    Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8); Nc = x(9); Na = x(10);
      
   
    dx = zeros(10,1);
    %dSc
    dx(1) = mu*Na - (betacc*Ic/Nc + betaac*Ia/Na + xicc*Ec/Nc +xiac*Ea/Na)*Sc - f*Sc;
    %dSa 
    dx(2) = f*Sc - (betaaa*Ia + betaca*Ic + xiaa*Ea + xica*Ec)*Sa - mu*Sa;
    %dEc
    dx(3) = (betacc*Ic/Nc + betaac*Ia/Na + xicc*Ec/Nc +xiac*Ea/Na)*Sc - epsilonc*Ec - f*Ec;
    %dEa
    dx(4) = (betaaa*Ia + betaca*Ic + xiaa*Ea + xica*Ec)*Sa - epsilona*Ea + f*Ec - mu*Ea;
    %dIc
    dx(5) = epsilonc*Ec - f*Ic - gammac*Ic;
    %dIa 
    dx(6) = f*Ic+ epsilona*Ea - gammaa*Ia - mu*Ia;
    %dRc
    dx(7) = gammac*Ic - f*Rc;
    %dRa
    dx(8) = gammaa*Ia +f*Rc - mu*Ra;
     %dNc
    dx(9) = mu*Na - f*Nc;
    %dNa
    dx(10) = f*Nc - mu*Na;

end


%%

function  error_in_data = err_in_dataSEIR_GLS_AdultChild_infectE(par, datac, dataa, T, noise)
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
Nc0 = Sc0 + Ic0 + Ec0 + Rc0;
Na0 = Sa0 + Ia0 + Ea0 + Ra0;
init_cond = [Sc0, Sa0, Ec0, Ea0, Ic0, Ia0, Rc0, Ra0, Nc0, Na0];

daily = 1:1:T; %Daily  

[t,y]=ode45(@SEIR_model_states_AdultChild_infectE, daily, init_cond,[],par); 

Infectedc=y(:,5);
Infecteda=y(:,6);
weightc=1./(noise* (1+Infectedc));  
weighta=1./(noise* (1+Infectedc)); 
error_in_data = sum(weightc.*(Infectedc-datac).^2) + sum(weighta.*(Infecteda-dataa).^2);
end