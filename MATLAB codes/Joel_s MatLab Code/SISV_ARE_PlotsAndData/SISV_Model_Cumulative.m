function [dx] = SISV_Model_Cumulative(t,x,z)
    betah = z(1);
    betav = z(2);
    gamma = z(3);
    muh = z(4);
    muv = z(5);
    pih = z(6);
    piv = z(7);


    Sh = x(1);Ih=x(2);Sv=x(3);Iv=x(4);Ch=x(5);Cv=x(6);
    Nh = Sh + Ih;
    Nv = Sv + Iv;
    dx = zeros(6,1);
    dx(1) = pih - muh*Sh-betah*Sh*Iv/Nv+gamma*Ih;
    dx(2) = betah*Sh*Iv/Nv-(muh + gamma)*Ih;
    dx(3) = piv - muv*Sv-betav*Sv*Ih/Nh;
    dx(4) = betav*Sv*Ih/Nh-muv*Iv;
    dx(5) = betah*Sh*Iv/Nv;
    dx(6) = betav*Sv*Ih/Nh;
end