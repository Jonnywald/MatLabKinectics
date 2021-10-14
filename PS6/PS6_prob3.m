%Problem Set 6
%prob 3

%RAW 6.6

%A -> B

% Species | Ci0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       |   3  |    Ca     |   3-z    |
% B       |   0  |    Cb     |    z     |
%---------|------|-----------|----------|
% Total   |   3  |    Ct     |    3     |

%r = k0 * exp(-E/T) * CA

E = 7550; %K
A = 3600; % cm2
Ta = 312; %K
T0 = 298; %K
Ca0 = 3.0; % Kmol/m3
U = 0.225; %J/cm2 s K
dHr = -8.09e8; %J/kmol
k0 = 4.48e6; %s-1
rho = 1e3; %Kg/m3
Cp = 4.19e3; %J/Kg/K
Vr = 18e-3; %m3
V0 = 60e-6; %m3/s

guess = [0.1 300];
result = fsolve(@(x)prob3(x,k0,E,U,A,Ta,Ca0,V0,Vr,T0,dHr,Cp,rho),guess);

Ca = result(1); %2.89 kmol/m3
T = result(2); %309.62 K

function f=prob3(x,k0,E,U,A,Ta,Ca0,V0,Vr,T0,dHr,Cp,rho)
    Ca = x(1);
    T = x(2);
    r = k0*exp(-E/T)*Ca;
    Q = U*A*(Ta - T);
    V = V0*(T/T0)*(Ca/Ca0);
    f(1) = Ca0*V0 - Ca*V + (-r)*Vr;
    f(2) = Q - rho*V0*Cp*(T0 - T) - dHr*(-r)*Vr;
end
