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

function f=prob2(x)

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


end
