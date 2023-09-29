function [dx]=compartmentsODEs(t,x,z) %time, iniital condition, parameters
%parameters
beta = z(1);
alpha = z(2);
k = z(3);
mod = z(4);

%inital condiitions


S = x(1); E = x(2); I = x(3); R = x(4);

    dx = zeros(4,1);
if mod == 0
    dx(1) = -beta*S*I;
    dx(3) = beta*S*I-alpha*I;
    dx(4) = alpha*I;
elseif mod == 1
    dx(1) = -beta*S*I;
    dx(2) = beta*S*I-k*E;
    dx(3) = k*E-alpha*I;
    dx(4) = alpha*I;
else
    print('no such model')
end

end


