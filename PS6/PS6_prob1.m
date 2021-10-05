%PROBLEM SET 6
%PROBLEM 1
%RAW 6.1

% A + B -> C

% Species | Ci0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       |   2  |    Ca     |  2 - z   |
% B       |   2  |    Cb     |  2 - z   |
% C       |   0  |    Cc     |  0 + z   |
%---------|------|-----------|----------|
% Total   |   4  |    Ct     |  4 - z   |

% r = kCaCb
% Ca = 2 - z
% z = 2 - Ca

% r = kCa^2

Ca0 = 2.0; % mol/L
Cb0 = 2.0; % mol/L
Cc0 = 0.0; % mol/L

Cpa = 20.0; % cal/mol.K
Cpb = 20.0; % cal/mol.K
Cpc = 40.0; % cal/mol.K

dHr = -10.0; % kcal/mol A @ 27C
VR = 1200; % L
Tr = 27.0; %C
Q = -41.7 * 1000; %cal/min
km = 0.01725; %L/mol/min
Tm = 300; %K
EaR = 2660; %K

%Unit conversions
Tr = Tr + 273.15; %converting C to K
dHr = dHr * 1000; %converting kcal/MolA to cal/MolA

%Material Balance
% dCadt = - kCa²

%Enrg. Balance

