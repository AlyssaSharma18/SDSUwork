function  error_in_data = err_in_dataSEIR_GLS(par,data)
%This is the error function for prevalence.    
%Initial Conditions
N = 1000;
I0 = 10;
E0 = 0;
R0 = 0;
S0 = N - I0 - E0 - R0;
init_cond = [S0,E0,I0,R0];



daily = 1:1:365; %Daily  


tspan=daily; 

[t,y]=ode45(@SEIR_model_states,daily,init_cond,[],par); 

Infected=y(tspan(:),3);
weight=1./Infected.^2;  
error_in_data = sum(weight.*(Infected-data).^2);

end