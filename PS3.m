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

%data given by the exercise:
k0 = 4.71e9; %s-1
Ea = 18.0; % Kcal/mol
V = 1000; % L
R = 1.98720425864083 / 1000; % Kcal * K-1 * mol-1
T = 60 + 273.15; % K
part1 = R*T;
part2 = -Ea/part1;
part3 = exp(part2);
k = k0 * part3; % s-1
k = k * 3600; % h-1

% densities of the components
densPO = 0.859;   % g/cm^3
densH2O = 1.000;  % g/cm^3
densPG = 1.0361;  % g/cm^3
densMet = 0.7914; % g/cm^3

%molar weigth of the components
MolWtPO = 58.08;  % g/mol
MolWtH2O = 18.02; % g/mol
MolWtPG = 76.11;  % g/mol
MolWtMet = 32.04; % g/mol

%volumetric feed rate of the components
VfPO = 1300;      % L/h
VfH2O = 6600;     % L/h
VfPG = 0;         % L/h
VfMet = 1300;     % L/h

%(a)
% for constant volume:
%total volumetric feed
VfT = VfPO + VfH2O + VfPG + VfMet;

%molar volume of the components
molVolPO = (MolWtPO / densPO)/1000; % mol / L
molVolH2O = (MolWtH2O / densH2O)/1000; % mol / L
molVolPG = (MolWtPG / densPG)/1000; % mol / L
molVolMet = (MolWtMet/densMet)/1000; % mol / L

%molar feed rate of each component
molFeedPO = VfPO * molVolPO; %Mol / h
molFeedH2O = molVolH2O * VfH2O; %Mol / h
molFeedPG = molVolPG * VfPG; %Mol / h
molFeedMet = molVolMet * VfMet; %Mol / h

%inital concentration of each component
cPO_0 = molFeedPO / VfPO; %mol / L
cH2O_0 = molFeedH2O / VfH2O; %mol / L
cPG_0 = 0;%molFeedPG / VfPG;(NaN) %mol / L
cMet_0 = molFeedMet / VfMet; %mol / L

%concentration of each component at the reactor exit (Ca)
cPO = cPO_0 / (1 + k * (V/VfT)); % 0.017480018491164
cH2O = cH2O_0 / (1 + k * (V/VfT)); % 0.004658683757370
cPG = cPO_0 - cPO; %0.050133485583341

%Conversion of PO
xPO = (cPO_0-cPO)/cPO_0; % 74.1471489602103%
%Production rate
rCPO = k*cPO; %0.461228067366740 mol / L * h

%for constant mass:

%total mass a the start
mass = (densMet * V / 0.001)/1000; %791.4 Kg

%(b)

%assuming all densities are equal to water
densH2O = 1.000;   % g/cm^3
densPO = densH2O;  % g/cm^3
densPG = densH2O;  % g/cm^3
densMet = densH2O; % g/cm^3

%molar volume of the components
molVolPO = (MolWtPO / densPO)/1000; % mol / L
molVolH2O = (MolWtH2O / densH2O)/1000; % mol / L
molVolPG = (MolWtPG / densPG)/1000; % mol / L
molVolMet = (MolWtMet/densMet)/1000; % mol / L

%molar feed rate of each component
molFeedPO = VfPO * molVolPO; %Mol / h
molFeedH2O = molVolH2O * VfH2O; %Mol / h
molFeedPG = molVolPG * VfPG; %Mol / h
molFeedMet = molVolMet * VfMet; %Mol / h

%inital concentration of each component
cPO_0 = molFeedPO / VfPO; %mol / L
cH2O_0 = molFeedH2O / VfH2O; %mol / L
cPG_0 = 0;%molFeedPG / VfPG;(NaN) %mol / L
cMet_0 = molFeedMet / VfMet; %mol / L

%concentration of each component at the reactor exit (Ca)
cPO = cPO_0 / (1 + k * (V/VfT)); % 0.015015335883910
cH2O = cH2O_0 / (1 + k * (V/VfT)); % 0.004658683757370
cPG = cPO_0 - cPO; %0.043064664116090

%Conversion of PO
xPO = (cPO_0-cPO)/cPO_0; % 0.741471489602103
%Production rate
rCPO = k*cPO; %0.396194909868030 mol / L * h

%smaller production rate when compared to using the "correct values"
%with a error of 14.1% in the production rate

%Problem 2
%plot concenration vs time of each component
%find time to max B

% A -> B <-> C

Ca0 = 1.0;
k1 = 0.5; %min-1
k2 = 0.25; %min-1
k3 = 0.1; %min-1

t_range = linspace(0,20,120);
global con kk1 kk2 kk3;
i = 1;
Ca = zeros(1,120);
Cb = zeros(1,120);
Cc = zeros(1,120);


for t = t_range
    Ca(i) = Ca0*exp(-k1*t);
    con = Ca(i);
    kk1 = k1;
    kk2 = k2;
    kk3 = k3;
    [y,x] = ode45(@prob2,t_range,[0 0]);
    Cb(i) = x(i,1);
    Cc(i) = x(i,2);
    i =i+1;
end


hold on
plot(t_range,Ca,"- k")
plot(t_range,Cb,"-- k")
plot(t_range,Cc,": k")
xlabel("time (min)")
ylabel("Concetration - Cj")
title("Concentration of components over time")
legend("Ca","Cb","Cc")
hold off

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

x0 = [0.9 1.0 1.1];
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
function u = prob2(t,x)
    global con kk1 kk2 kk3;
    dcbdt = kk1*con - kk2*x(1) + kk3*x(2);
    dccdt = kk2*x(1) - kk3*x(2);
    u = [dcbdt;dccdt];
end