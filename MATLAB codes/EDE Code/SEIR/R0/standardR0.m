function y = standardR0(par)
    beta = par(1);
    epsilon = par(2);
    gamma = par(3);
    mu = 0.00007;
    N = 1000;

    y = beta*epsilon*N/((gamma + mu)*(epsilon + mu));


end
