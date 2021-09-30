function rdxdt=myfunc(t,x)
    k1 = 0.02;
    k2 = 0.00002;
    k3 = 0.00004;
    k4 = 0.04;
    dxdt = k1*x(1) - k2*x(1)*x(2);
    dydt = k3*x(1)*x(2) - k4*x(2);
    rdxdt = [dxdt; dydt];
end