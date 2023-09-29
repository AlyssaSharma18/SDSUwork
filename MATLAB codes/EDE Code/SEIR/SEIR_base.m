T=365;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.00027 0.00027*250 0.00027*750 0.00027 0.087 0.087 0.25 0.25]; %betaaa, betacc, betaac, betaca, gammac, gammaa, epsilonc, epsilona
Fitted_Parameters2 = [0.00027 0.25 0.087]; %beta, epsilon, gamma
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


S0 = 980;
E0 = 0;
I0 = 20;
R0 = 0;
N0 = S0 + E0 + I0 + R0;
init_cond2= [S0,E0,I0,R0];


% Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states_AdultChild_a,daily,init_cond,[],Fitted_Parameters);

figure; hold on
plot(t,y_est(:,1)+y_est(:,2))
plot(t,y_est(:,3)+y_est(:,4))
plot(t,y_est(:,5)+y_est(:,6))
plot(t,y_est(:,7)+y_est(:,8))
yline(50)
yline(100)
yline(150)
yline(200)
yline(250)
yline(500)
yline(600)
yline(700)
yline(750)
yline(800)

hold off

[t,y_est2] = ode45(@SEIR_model_states,daily,init_cond2,[],Fitted_Parameters2);

figure
hold on
plot(t,y_est2)

yline(50)
yline(100)
yline(150)
yline(200)
yline(250)
yline(500)
yline(600)
yline(700)
yline(750)
yline(800)

hold off
%%  SEIR Adult Child model
function [dx]  = SEIR_model_states_AdultChild_a(t,x,z)

    %parameters to be estimated
    %beta = z(1); %0.00027
    betaaa = z(1);
    betacc = z(2); %(1/4)*(24515/42739)*250*beta;
    betaac = z(3); %(1/6)*750*beta;
    betaca = z(4); %(1/3)*(24515/42739)*beta;
    gammaa = z(5); %0.074; 
    gammac = z(6); %0.1; 
    epsilonc = z(7); %0.3
    epsilona = z(8); %0.2
    mu = 0; %0.00007;
    f = 0; %0.0005;
    Nc0 = 250;
    Na0 = 750;
      
    Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8); Nc = x(9); Na = x(10);
    
   %IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0
   
   
    dx = zeros(10,1);
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
    %dNc
    dx(9) = mu*Na - f*Nc;
    %dNa
    dx(10) = f*Nc - mu*Na;


    
end

%% SEIR Standard Model
function [dx]  = SEIR_model_states(t,x,z)

    %parameters to be estimated
    beta = z(1);
    epsilon = z(2);
    gamma = z(3);

      
    S = x(1); E = x(2); I = x(3); R= x(4);
    
   %IC: Total N=1000. S(0)=990 E(0)=0 I(0)=10  R(0)=0
   
   
    dx = zeros(4,1);
    %dS
    dx(1) = -beta*S*I;
    %dE 
    dx(2) = beta*S*I-epsilon*E;
    %dI
    dx(3) = epsilon*E-gamma*I;
    %dR
    dx(4) = gamma*I;

    
end

