%% Calculating Ro for the SEIR Adult-Child model
syms xi_cc xi_ac xi_ca xi_aa beta_cc beta_ac beta_ca beta_aa S_c S_a ...
    epsilon_c epsilon_a mu f gamma_c gamma_a N_c N_a
% F = [xi_cc*S_c/N_c, xi_ac*S_c/N_a, beta_cc*S_c/N_c, beta_ac*S_c/N_a; ...
% xi_ca*S_a, xi_aa*S_a, beta_ca*S_a, beta_aa*S_a; 0, 0, 0, 0; 0, 0, 0, 0];

V = [epsilon_c + f, 0, 0, 0; -f, epsilon_a + mu, 0, 0; -epsilon_c, 0, f + gamma_c, 0;...
    0, -epsilon_a, -f, gamma_a + mu];

W = inv(V);

newF = [0, 0,beta_cc, beta_ac*N_c/N_a; 0, 0, beta_ca*N_a, beta_aa*N_a;... 
    0, 0, 0, 0; 0, 0, 0, 0];
newNGM = newF*W;
newprange = eig(newNGM);

% dbaa = diff(newprange,beta_aa);

% NGM = F*W;
% prange = eig(NGM);

 syms a b c d e f g h

A = beta_cc*epsilon_c/((epsilon_c+f)*(f+gamma_c)) + N_c*beta_ac*((epsilon_a*f^2 ...
    + epsilon_a*epsilon_c*f + epsilon_a*f*gamma_c + epsilon_c*f*mu)/(N_a*(epsilon_c+f)*(f+ ...
    gamma_c)*(epsilon_a+mu)*(gamma_a+mu)));
B = N_c*beta_ac*(epsilon_a/(N_a*(epsilon_a + mu)*(gamma_a + mu)));
E = N_a*(2*beta_ca)*(epsilon_c/((epsilon_c+f)*(f+gamma_c))) + N_a*beta_aa*((epsilon_a*f^2 ...
    + epsilon_a*epsilon_c*f + epsilon_a*f*gamma_c + epsilon_c*f*mu)/((epsilon_c ...
    + f)*(f + gamma_c)*(epsilon_a + mu)*(gamma_a + mu)));
F = N_a*beta_aa*(epsilon_a/((epsilon_a+mu)*(gamma_a+mu)));

%  Sample = [a, b, c, d; e, f, g, h; 0, 0, 0, 0; 0, 0, 0, 0];
%  samplerange = eig(Sample);

Ro = A/2 + F/2 + ((A^2 + F^2 - 2*A*F + 4*B*E)^(1/2))/2;

depsilon_c = diff(Ro,epsilon_c);

dec = double(subs(depsilon_c, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
                                {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.5 0.5 0.000000005 0.00003}));

depsilon_a = diff(Ro,epsilon_a);

dea = double(subs(depsilon_a, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.5 0.5 0.000000005 0.00003}));

dgamma_a= diff(Ro,gamma_a);

dga = double(subs(dgamma_a, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));

dgamma_c = diff(Ro,gamma_c);

dgc = double(subs(dgamma_c, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));

dbeta_aa = diff(Ro,beta_aa);

db_aa_value = double(subs(dbeta_aa, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));

dbeta_ac = diff(Ro,beta_ac);

db_ac_value = double(subs(dbeta_ac, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));

dbeta_cc = diff(Ro,beta_cc);

db_cc_value = double(subs(dbeta_cc, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));

dbeta_ca = diff(Ro,beta_ca);

db_ca_value = double(subs(dbeta_ca, {beta_aa beta_ac beta_cc beta_ca N_c N_a gamma_a gamma_c epsilon_c epsilon_a f mu}, ...
    {0.0019 0.013 0.18 0.0001 250 750 0.0003 0.04 0.05 0.05 0.000000005 0.00003}));



