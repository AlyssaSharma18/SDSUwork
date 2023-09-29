% Date Created: 7/12/2022
% Created by Anuradha Agarwal
%This is to test the gamma with non dependent beta and xi 
%Time span
T=200;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.00027 0.01 0.01 0.00005 0.1 0.074 0.3 0.2 0.0067 0.00018 0.000033 0.0067]; 

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

[t,y_est]=ode45(@SIR_model_states,daily,init_cond,[],Fitted_Parameters);


Prevalencedatac = y_est(:,5); %+ 0.00*(y_est(:,5)) ;
Prevalencedataa = y_est(:,6); %0.00*(y_est(:,6));
% figure
% hold on
% plot(t,y_est(:,5));
% plot(t,y_est(:,6))
% plot(t,Prevalencedatac)
% plot(t,Prevalencedataa)


% [t,y_esta]=ode45(@SIR_model_states,daily,init_cond,[],paramtest);
figure
hold on
plot(t,y_est(:,5), 'blue')
plot(t,y_est(:,6), 'red')
legend('Ic','Ia')
hold off
% plot(t,y_est(:,6))
% plot(t,y_esta(:,5))
% plot(t,y_esta(:,6))

%% MC Optimization
parameter1 = linspace(0, 0.1, 50);
parameter2 = linspace(0.11, 0.3, 50);
firstparameter = [parameter1 parameter2];
length_of_beta = length(firstparameter);
EstiParam = zeros(length_of_beta, Number_Parameters); %Storage for estimated parameters.
function_values = zeros(1, length_of_beta);
beta_val = zeros(1,length_of_beta);
alpha_val = zeros(1,length_of_beta);
j = [0.00027 0.01 0.01 0.00005 0.074 0.3 0.2 0.0067 0.00018 0.000033 0.0067];
%this is where you change the parameter (delete it)
for i = 1:length_of_beta
    Initial_Guess = j;
    
    %Bounds for the parameters
    Lowerbounds = [0 0 0 0 0 0 0 0 0 0 0];
    Upperbounds = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    fixp = firstparameter(i);
    options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
    [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_data_SIR(Initial_Guess,Prevalencedatac, Prevalencedataa, T,fixp),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);
    EstiParam(i,:) = [fixp,EstimatedParameters]; 
    function_values(1,i) = fval;
    % j = EstimatedParameters;
end

%%

thresholdup = chi2inv(0.95,12) + min(function_values);
threshold = chi2inv(0.95,1) + min(function_values);
figure
plot(firstparameter, function_values)
hold on 
yline(threshold, 'r','Threshold')
yline(thresholdup, 'r', 'Threshold max' )
title('Profile Likelihood {\gamma_{C}}')
hold off 


%%  SEIR model
function [dx]  = SIR_model_states(t,x,z)


betaaa = z(1); %5*10^-7 to 2*10^-
betacc = z(2); %0.1 to 1
betaac = z(3); %.05 to 1 <betacc
betaca = z(4); % 1*10^-8 to 5*10^-7
gammac = z(5); %.02 to 1
gammaa = z(6); % .01 to 1
epsilonc = z(7); % .04 to 1
epsilona = z(8); % .02 to 1
xicc = z(9); %0.2
xiaa = z(10); %0.000001
xica = z(11); %0.1
xiac = z(12);%0.00000002;
mu = 0.00008; %4*10^-5 to 8*10^-5
f = 0.0005; %1*10^-4 to 5*10^-4
Nc = 250;
Na = 750;


      
Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8);
    
   %IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0
   
   
     dx = zeros(8,1);
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
%     %dNc
%     dx(9) = mu*Na - f*Nc;
%     %dNa
%     dx(10) = f*Nc - mu*Na;

end

%%

function  error_in_data = err_in_data_SIR(par,datac,dataa,T,fp)
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

daily = 1:1:T; %Daily
tspan=daily;
par1=[par(1:4), fp, par(5:end)];;
[t,y]=ode45(@SIR_model_states,daily,init_cond,[],par1);


Infectedc=y(tspan(:),5);
Infecteda=y(tspan(:),6);
weightc=1;  
weighta=1; 
error_in_data = sum((weightc.*(Infectedc-datac)).^2) + sum((weighta.*(Infecteda-dataa)).^2);

end



