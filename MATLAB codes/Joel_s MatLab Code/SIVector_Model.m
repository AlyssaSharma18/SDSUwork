function [dx] = SIVector_Model(t,x,z)
    muh = z(1);
    muv = z(2);
    betah = z(3);
    betav = z(4);
    gamma = z(5);
    pih = z(6);
    piv = z(7);


    Sh = x(1);Ih=x(2);Sv=x(3);Iv=x(4);

    dx = zeros(4,1);
    dx(1) = pih - muh*Sh-betah*Sh*Iv+gamma*Ih;
    dx(2) = betah*Sh*Iv-(muh + gamma)*Ih;
    dx(3) = piv - muv*Sv-betav*Sv*Ih;
    dx(4) = betav*Sv*Ih-muv*Iv;
end