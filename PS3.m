%PS3

% Problem set 3

% Guilherme Bertola

% Problem 1
% RAW 4.8

% propilene Oxide + H2O -> propylene glycol

%to simplify: PO + H2O -> PG

%excess of H2O

% r = k * CPO

% k = k0 * exp((-Ea)/R*T)
k0 = 4.71e9; %s-1
Ea = 18.0; % Kcal/mol
V = 1000; % L
R = 1.98720425864083 / 1000; % Kcal * K-1 * mol-1
T = 60 + 273.15; % K
k = k0 * exp((-Ea)/R*T); % s-1

densPO = 0.859;   % g/cm^3
densH2O = 1.000;  % g/cm^3
densPG = 1.0361;  % g/cm^3
densMet = 0.7914; % g/cm^3

MolWtPO = 58.08;  % g/mol
MolWtH20 = 18.02; % g/mol
MolWtPG = 76.11;  % g/mol
MolWtMet = 32.04; % g/mol

VfPO = 1300;      % L/h
VfH2O = 6600;     % L/h
VfPG = 0;         % L/h
VfMet = 1300;     % L/h


