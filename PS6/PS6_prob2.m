%problem set 6
%problem 2

%RAW 6.5

% A -> B (1)
% A -> C (2)

% Species | Ci0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       |   5  |    Ca     | 5-z1-z2  |
% B       |   0  |    Cb     |    z1    |
% C       |   0  |    Cc     |    z2    |
%---------|------|-----------|----------|
% Total   |   5  |    Ct     |    5     |

Ca0 = 5.0; %mol / L
Vr = 1000; % L
V0 = 100; % L/min
Cp = 0.22; %cal/(g K)
dHr1 = -12e3; %Cal/Mol
dHr2 = -15e3; %cal/Mol
Mf = 93.2e3; %g/min

T0 = 45 + 273.15; %K

%a)
Selectivity = 5;

guess = [1.0 2.0 0.1 400 -3000];
result = fsolve(@(x)prob2(x,Mf,Cp,Ca0,V0,Vr,dHr1,dHr2,Selectivity,T0),guess);

Q = result(5); % -3.821102233701143e+06 ??????

%b)
Selectivity = 4;
guess = [1.0 2.0 0.1 400 -3000];
result = fsolve(@(x)prob2(x,Mf,Cp,Ca0,V0,Vr,dHr1,dHr2,Selectivity,T0),guess);

Q = result(5); % -1.130489902080836e+06 ??????

%c)
%in order to achive a greater selectivity it is needed to remove a larger
%amont of energery, that can be percieved in the Q being a larger negative
%value when the seletivity desired is 5 instead of 4

function f=prob2(x,Mf,Cp,Ca0,V0,Vr,dHr1,dHr2,Selectivity,T0)

Ca = x(1);
Cb = x(2);
Cc = x(3);
T = x(4);
Q = x(5);

k1 = 3.16e14 * exp(-12500/T);
k2 = 2.52e9 * exp(-8500/T);

r1 = k1 * Ca;
r2 = k2 * Ca;

ra = -r1 -r2;
rb = r1;
rc = r2;


V = V0*(T/T0)*(Ca/Ca0);

f(1) = Ca0*V0 - Ca*V + ra*Vr;
f(2) = -Cb*V + rb*Vr;
f(3) = -Cc*V + rc*Vr;
f(4) = Mf*Cp*(T0 - T) + Q - (dHr1*r1 + dHr2*r2)*Vr;
f(5) = rb - Selectivity*rc;

end
