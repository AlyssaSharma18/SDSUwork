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
