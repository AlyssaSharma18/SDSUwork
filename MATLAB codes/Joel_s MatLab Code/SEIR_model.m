function [dx] = SEIR_model(t,x,z)
    beta = z(1);
    gamma = z(2);
    alpha = z(3);

    S = x(1);E=x(2); I=x(3);R=x(4);

    dx = zeros(4,1);

    dx(1) = -beta*S*I;
    dx(2) = beta*S*I-gamma*E;
    dx(3) = gamma*E-alpha*I;
    dx(4) = alpha*I;
end