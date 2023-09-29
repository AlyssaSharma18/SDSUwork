function [dx] = SEIR_model_states_AdultChild(t,x,z)

%parameters to be estimated

betaaa = z(1);
betacc = z(2);
betaac = z(3);
betaca = z(4);
gammac = z(5);
gammaa = z(6);
epsilonc = z(7);
epsilona = z(8);
mu = 0.00007;
f = 0.00005;
Nc = 250;
Na = 750;


Sc = x(1); Sa = x(2); Ec = x(3); Ea = x(4); Ic = x(5); Ia = x(6); Rc = x(7); Ra = x(8);

dx = zeros(8,1);
%dSc
dx(1) = mu*Na - (betacc*Ic/Nc + betaac*Ia/Na)*Sc - f*Sc;
%dSa
dx(2) = f*Sc - (betaaa*Ia + betaca*Ic)*Sa - mu*Sa;
%dEc
dx(3) = (betacc*Ic/Nc + betaac*Ia/Na)*Sc - epsilonc*Ec - f*Ec;
%dEa
dx(4) = (betaaa*Ia + betaca*Ic)*Sa - epsilona*Ea + f*Ec - mu*Ea;
%dIc 
dx(5) = epsilonc*Ec - f*Ic - gammac*Ic;
%dIa 
dx(6) = f*Ic + epsilona*Ea - gammaa*Ia - mu*Ia;
%dRc
dx(7) = gammac*Ic - f*Rc;
%dRa
dx(8) = gammaa*Ia + f*Rc - mu*Ra;

end

