function [dx] = SIR_model(t,x,z)
    beta = z(2);
    alpha = z(1);

    S = x(1);I=x(2);R=x(3);

    dx=zeros(3,1);
    dx(1) = -beta*S*I;
    dx(2) = beta*S*I-alpha*I;
    dx(3) = alpha*I;
end