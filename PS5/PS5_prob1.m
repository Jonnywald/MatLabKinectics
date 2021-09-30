%Problem set 5
%problem 1

%RAW 6.2

% A + B -> C

% r = kCaCb

Ca0 = 2.0; % mol/L
Cb0 = 2.0; % mol/L
Cc0 = 0.0; % mol/L

km = 0.01725; % L/mol.min
dHr = -10.0; % kcal/mol A @ 27C
Cpa = 20.0; % cal/mol.K
Cpb = 20.0; % cal/mol.K
Cpc = 40.0; % cal/mol.K

Tm = 300.0; %K
EaR = 2660.0; %K

VR = 1200; % L

T1 = 27.0; %C

Q = -41.7 * 1000; %cal/min

% Enrg balance
% 0 = Q - (Fa0*Cpa(T-T0)+Fb0*Cb0(T-T0)) - (dHr*dCp(T-Tr))

%Unit conversions
T1 = T1 + 273.15; %converting C to K
dHr = dHr * 1000; %converting kcal/MolA to cal/MolA


