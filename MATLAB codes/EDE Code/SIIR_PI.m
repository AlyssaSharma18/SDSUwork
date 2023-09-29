clear all
%Number of iterations for the simulation.
NumberofIterations=2;

%%
%Time span
T=365;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.02 0.01 0.05 0.05 0.04 0.02];%[0.0000005 0.1 0.05 0.00000001 0.02 0.01 0.04 0.02 0.2 0.000001 0.1 0.00000002]; %betaaa, betacc, betaac, betaca, gammac, gammaa, epsilonc, epsilona
%Peak at day 109.

%Initial Conditions
Sc0 = 240;
Sa0 = 740;
Iqc0 = 10;
Iqa0 = 10;
Ic0 = 0;
Ia0 = 0;
Rc0 = 0;
Ra0 = 0;
Nc0 = Sc0 + Ic0 + Iqc0 + Rc0;
Na0 = Sa0 + Ia0 + Iqa0 + Ra0;
init_cond = [Sc0, Sa0, Ic0, Ia0, Iqc0, Iqa0, Rc0, Ra0, Nc0, Na0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters.

ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.
%Fval = zeros(1,NumberofIterations); %Storage for optimization values.

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters.
Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels from paper.

%Uses estimated parameters to solve system.
[t,y_est]=ode45(@SIIR_model_states_AdultChild_infectE,daily,init_cond,[],Fitted_Parameters);


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
        Lowerbounds = [0 0 0 0 0 0];
        Upperbounds=[1 1 1 1 1 1];

        %Optimization
        %Fminsearchbnd
        %[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Initial_Guess)err_in_dataSEIR_GLS(Initial_Guess,Prevalencedata),Initial_Guess, Lowerbounds, Upperbounds, options);
        options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
        [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_dataSIIR_GLS_AdultChild_infectE(Initial_Guess,Prevalencedatac,Prevalencedataa, T),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);


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
function [dx]  = SIIR_model_states_AdultChild_infectE(t,x,z)

    %parameters to be estimated
    betaaa = 0.0000005; %5*10^-7 to 2*10^-5
    betacc = 0.1; %0.1 to 1
    betaac = 0.05; %.05 to 1 <betacc
    betaca = 0.00000001; % 1*10^-8 to 5*10^-7
    gammaic = z(1); %.02 to 1
    gammaia = z(2); % .01 to 1
    gammaqc = z(3); %0.05;
    gammaqa = z(4); %0.05;
    ic = z(5); % .04 to 1
    ia = z(6); % .02 to 1
    xicc = 0.2; %0.2
    xiaa = 0.000001; %0.000001
    xica = 0.1; %0.1
    xiac = 0.00000002;%0.00000002;
    mu = 0.00008; %4*10^-5 to 8*10^-5
    f = 0.0005; %1*10^-4 to 5*10^-4
    Nc0 = 250;
    Na0 = 750;


      
    Sc = x(1); Sa = x(2); Ic = x(3); Ia = x(4); Iqc = x(5); Iqa = x(6); Rc = x(7); Ra = x(8); Nc = x(9); Na = x(10);
    
   %IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0
   
   
    dx = zeros(10,1);
    %dSc
    dx(1) = mu*Na - (betacc*Iqc/Nc + betaac*Iqa/Na + xicc*Ic/Nc +xiac*Ia/Na)*Sc - f*Sc;
    %dSa 
    dx(2) = f*Sc - (betaaa*Iqa + betaca*Iqc + xiaa*Ia + xica*Ic)*Sa - mu*Sa;
    %dEc
    dx(3) = (betacc*Iqc/Nc + betaac*Iqa/Na + xicc*Ic/Nc +xiac*Ia/Na)*Sc - ic*Ic - gammaic*Ic - f*Ic;
    %dEa
    dx(4) = (betaaa*Iqa + betaca*Iqc + xiaa*Ia + xica*Ic)*Sa - ia*Ia + f*Ic - mu*Ia - gammaia*Ia;
    %dIc
    dx(5) = ic*Ic - f*Iqc - gammaqc*Iqc;
    %dIa 
    dx(6) = f*Iqc+ ia*Ia - gammaqa*Iqa - mu*Iqa;
    %dRc
    dx(7) = gammaqc*Iqc - f*Rc;
    %dRa
    dx(8) = gammaqa*Iqa +f*Rc - mu*Ra;
     %dNc
    dx(9) = mu*Na - f*Nc;
    %dNa
    dx(10) = f*Nc - mu*Na;

    
end


%%

function  error_in_data = err_in_dataSIIR_GLS_AdultChild_infectE(par, datac, dataa, T)
%This is the error function for prevalence.    
%Initial Conditions
Sc0 = 240;
Sa0 = 740;
Iqc0 = 10;
Iqa0 = 10;
Ic0 = 0;
Ia0 = 0;
Rc0 = 0;
Ra0 = 0;
Nc0 = Sc0 + Ic0 + Iqc0 + Rc0;
Na0 = Sa0 + Ia0 + Iqa0 + Ra0;
init_cond = [Sc0, Sa0, Ic0, Ia0, Iqc0, Iqa0, Rc0, Ra0, Nc0, Na0];



daily = 1:1:365; %Daily  


tspan=daily; 

[t,y]=ode45(@SIIR_model_states_AdultChild_infectE, daily, init_cond,[],par); 

Infectedc=y(tspan(:),5);
Infecteda=y(tspan(:),6);
weightc=1./Infectedc.^2;  
weighta=1./Infecteda.^2; 
error_in_data = sum(weightc.*(Infectedc-datac).^2) + sum(weighta.*(Infecteda-dataa).^2);
end