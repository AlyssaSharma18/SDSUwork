T=120;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters= [0.00027 0.074 0.1 0.3 0.2]; %betaaa, betacc, betaac, betaca, gammaa, gammac, epsilonc, epsilona
Fitted_Parameters2 = [0.00040536 0.25 0.087];%[0.00005 0.3 0.1];%[0.000253 0.25 0.087]; %beta, epsilon, gamma
%Peak at day 109.

%Initial Conditions
Sc0 = 250;
Sa0 = 750;
Ic0 = 10;
Ia0 = 10;
Ec0 = 0;
Ea0 = 0;
Rc0 = 0;
Ra0 = 0;
Nc0 = Sc0 + Ic0 + Ec0 + Rc0;
Na0 = Sa0 + Ia0 + Ea0 + Ra0;
init_cond = [Sc0, Sa0, Ec0, Ea0, Ic0, Ia0, Rc0, Ra0, Nc0, Na0];


S0 = 1000;%1000;
E0 = 0;
I0 = 10; %20
R0 = 0;
N0 = S0 + E0 + I0 + R0;
init_cond2= [S0,E0,I0,R0];


% Uses estimated parameters to solve system.
[t,y_est]=ode45(@SEIR_model_states_AdultChild_a,daily,init_cond,[],Fitted_Parameters);
[t,y_est2] = ode45(@SEIR_model_states,daily,init_cond2,[],Fitted_Parameters2);

figure; hold on
% plot(t,y_est(:,1))
% plot(t,y_est(:,2))
% plot(t,y_est(:,3))
% plot(t,y_est(:,4))
plot(t,y_est(:,5))
plot(t,y_est(:,6))
% plot(t,y_est(:,7))
% plot(t,y_est(:,8))
% plot(t,y_est2(:,3), '--')
hold off
%legend('I_A','I')
%legend('S_A','E_A','I_A','R_A')
%saveas(gcf,'SEIRsolutions.jpg')
%legend('S_C','S_A','E_C','E_A','I_C','I_A','R_C','R_A')
% title('Adult Child R_0 = 2.0128, Adult SEIR R_0 = 2.7401')
%saveas(gcf,'AdultSEIR_infectiousChild.jpg')

% figure; hold on
% plot(t,y_est(:,1))
% plot(t,y_est(:,3))
% plot(t,y_est(:,5))
% plot(t,y_est(:,7))
% hold off
% legend('S_C','E_C','I_C','R_C')
%saveas(gcf,'ChildSolutions.jpg')
% 
% figure; hold on
% plot(t,y_est(:,2))
% plot(t,y_est(:,4))
% plot(t,y_est(:,6))
% plot(t,y_est(:,8))
% hold off
% legend('S_A','E_A','I_A','R_A')
% saveas(gcf,'AdultSolutions.jpg')

% figure; hold on
% plot(t,y_est(:,1)+y_est(:,2))
% plot(t,y_est(:,3)+y_est(:,4))
% plot(t,y_est(:,5)+y_est(:,6))
% plot(t,y_est(:,7)+y_est(:,8))
% 
% 
% plot(t,y_est2,'--')
% 
% hold off
% % % 
% legend('S_C + S_A','E_C + E_A', 'I_C + I_A', 'R_C + R_A', 'S', 'E', 'I', 'R')
% title('SEIR R_0 = 4.6593 ; Adult Child R_{0} = 2.0128')
% 
%saveas(gcf,'infectiousChildmidepidemic.jpg')
%%  SEIR Adult Child model
function [dx]  = SEIR_model_states_AdultChild_a(t,x,z)

    %parameters to be estimated
    beta = z(1); %0.00027
    betaaa = beta; %beta/750; %z(1);
    betacc = (1/4)*(24515/42739)*250*beta;
    betaac = (1/6)*750*beta;
    betaca = (1/3)*(24515/42739)*beta;
    gammaa = z(2); %0.074; 
    gammac = z(3); %0.1; 
    epsilonc = z(4); %0.3
    epsilona = z(5); %0.2
    mu = 0.00007;
    f = 0.0005;
    Nc0 = 260;
    Na0 = 760;
      
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

