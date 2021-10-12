%problem set 6
%problem 2

%RAW 6.5

% A -> B (1)
% A -> C (2)

Ca0 = 5.0; %mol / L
Vr = 1000; % L
V0 = 100; % L/min
Cp = 0.22; %cal/(g K)
dHr1 = -12e3; %Cal/Mol
dHr2 = -15e3; %cal/Mol
Mf = 93.2e3; %g/min

T0 = 45 + 273.15; %K

Selectivity = 5;

guess = [1.0 2.0 2.0 400];
result = fsolve(@(x)prob2(x,Mf,Cp,Ca0,V0,Vr,dHr1,dHr2,Selectivity,Ta,T0),guess);

function f=prob2(x,Mf,Cp,Ca0,V0,Vr,dHr1,dHr2,Selectivity,Ta,T0)

Ca = x(1);
Cb = x(2);
Cc = x(3);
T = x(4);


k1 = 3.16e14 * exp(-12500/T);
k2 = 2.52e9 * exp(-8500/T);

r1 = k1 * Ca;
r2 = k2 * Ca;

ra = -r1 -r2;
rb = r1;
rc = r2;

Q = Mf*Cp*(Ta - T);

V = V0*(T/T0)*(Ca/Ca0);

f(1) = Ca0*V0 - Ca*V + ra*Vr;
f(2) = -Cb*V + rb*Vr;
f(3) = -Cc*V + rc*Vr;
f(4) = Q -Ca0*V*Cp*(T - T0) - (dHr1 + dHr2)*ra*Vr;
f(5) = rb - Selectivity*rc;

end
