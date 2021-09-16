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

%Problem 4
% use PREOS to calculate the fulgacity of CO2 and CH4 in a 50/50 Mixute at 
%500K and 500 bar

T = 500; %K
Pt = 500; %Bar
R = 0.0831446261815324; %L*ber*k-1*mol-1

MolWtCO2 = 44.010; %g/mol
TcCO2 = 304.2; %K
PcCO2 = 7.376 * 10; %bar
OmegaCO2 = 0.225; %accentric factor
yCO2 = 0.5; % molar fraction

MolWtCH4 = 16.04; %g/mol
TcCH4 = 190.55; %K
PcCH4 = 4.595 * 10; %bar
OmegaCH4 = 0.01142; %accentric factor
yCH4 = 0.5; % molar fraction

kCO2 = 0.37464 +1.54226*OmegaCO2 - 0.26992*(OmegaCO2.^2);
kCH4 = 0.37464 +1.54226*OmegaCH4 - 0.26992*(OmegaCH4.^2);

alfaCO2 = 1 + kCO2*(1-sqrt(T/TcCO2));
alfaCH4 = 1 + kCH4*(1-sqrt(T/TcCH4));

aCO2 = 0.45724*(((R.^2)*(TcCO2.^2))/PcCO2)*alfaCO2;
aCH4 = 0.45724*(((R.^2)*(TcCH4.^2))/PcCH4)*alfaCH4;

bCO2 = 0.07780*(R*TcCO2)/(PcCO2);
bCH4 = 0.07780*(R*TcCH4)/(PcCH4);

ACO2 = (aCO2*Pt)/((R*T).^2);
ACH4 = (aCH4*Pt)/((R*T).^2);

BCO2 = (bCO2*Pt)/(R*T);
BCH4 = (bCH4*Pt)/(R*T);

aCO2CH4 = sqrt(aCO2*aCH4);
ACO2CH4 = (aCO2CH4*Pt)/((R*T).^2);

aMix = yCO2*yCO2*aCO2 + yCH4*yCH4*aCH4 + 2*(yCO2*yCH4*aCO2CH4);
bMix = yCO2*bCO2 + yCH4*bCH4;

AMix = (aMix*Pt)/((R*T).^2);
BMix = (bMix*Pt)/(R*T);

x0 = [0.09 0.2 0.9];
z = fsolve(@(z)cmp(z,AMix,BMix),x0);
Z = max(z);

part1_1 = (BCO2/BMix)*(Z-1);
part1_2 = (BCH4/BMix)*(Z-1);

part2_1 = log(Z-BMix);
part2_2 = log(Z-BMix);

part3_1 = (AMix/(2*sqrt(2)*BMix));
part3_2 = (AMix/(2*sqrt(2)*BMix));

part4_1 = (((2*(yCO2*ACO2 + yCH4*ACO2CH4))/AMix)-((BCO2/BMix)));
part4_2 = (((2*(yCH4*ACH4 + yCO2*ACO2CH4))/AMix)-((BCH4/BMix)));

part5_1 = log((Z+(1+sqrt(2))*BMix)/(Z+(1-sqrt(2))*BMix));
part5_2 = log((Z+(1+sqrt(2))*BMix)/(Z+(1-sqrt(2))*BMix));

lnphiCO2 = part1_1 - part2_1 - part3_1*part4_1*part5_1;
lnphiCH4 = part1_2 - part2_2 - part3_2*part4_2*part5_2;

phiCO2 = exp(lnphiCO2);
phiCH4 = exp(lnphiCH4);

fulgacityCO2 = phiCO2 * yCO2 * Pt; %f = 180.1616463507124 bar
fulgacityCH4 = phiCH4 * yCH4 * Pt; %f = 241.9417151382955 bar

function f = cmp(z,A,B)
    f = z.^3 - (1-B)*z.^2 + (A-3*(B.^2)-2*B)*z - (A*B -B.^2 - B.^3);
end
