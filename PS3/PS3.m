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

% Species  | Feed |    Exit (Ci)    |               yi               |     
%--------------------------------------------------------------------|
% PO       | Cpo0 |    Cpo0 - z     | (Cpo0 - z)/(Cpo0+Ch2o0+Cmet-z) |
% H2O      | Ch2o0|    Ch2o0 - z    | (Ch2o0 - z)/(Cpo0+Ch2o0+Cmet-z)|
% PG       |  0   |        z        |    (z)/(Cpo0+Ch2o0+Cmet-z)     |
% Met      | CMet |      Cmet       |    Cmet/(Cpo0+Ch2o0+Cmet-z)    |
%----------|------|-----------------|--------------------------------|
% Total    |  Ct  |Cpo0+Ch2o0+Cmet-z|                1               |

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
molVolPO = (MolWtPO / densPO)/1000; % L / mol
molVolH2O = (MolWtH2O / densH2O)/1000; % L / mol
molVolPG = (MolWtPG / densPG)/1000; % L / mol
molVolMet = (MolWtMet/densMet)/1000; % L / mol

%molar feed rate of each component
molFeedPO = VfPO / molVolPO; %Mol / h
molFeedH2O = VfH2O / molVolH2O; %Mol / h
molFeedPG = VfPG / molVolPG; %Mol / h
molFeedMet = VfMet / molVolMet; %Mol / h

%inital concentration of each component
cPO_0 = molFeedPO / VfPO; %mol / L
cH2O_0 = molFeedH2O / VfH2O; %mol / L
cPG_0 = 0;%molFeedPG / VfPG;(NaN) %mol / L
cMet_0 = molFeedMet / VfMet; %mol / L

%concentration of each component at the reactor exit (Ca)
cPO = cPO_0 / (1 + k * (V/VfT)); % 3.823622424789828 Mol / L
cH2O = cH2O_0 / (1 + k * (V/VfT)); % 14.346754184123006 Mol / L
cPG = cPO_0 - cPO; %10.966322478791440 Mol / L

%Conversion of PO
xPO = (cPO_0-cPO)/cPO_0; % 74.1471489602103%
%Production rate
rCPO = k*cPO; %100.8901668048813 mol / L * h

%for constant mass:

%total mass at the start
mass = (densMet * V / 0.001); %791400 g

% V = Mass/density
% V = Mass/(Cpo*Mpo+Ch2o*Mh2o+Cpg*Mpg+Cmet*Mmet)
% 0 = Mass/(Cpo*Mpo+Ch2o*Mh2o+Cpg*Mpg+Cmet*Mmet) - V (1)

% 1 = Cpo*V0po + Ch2o*V0h2o + Cpg*V0pg + Cmet*V0met
% 0 = Cpo*V0po + Ch2o*V0h2o + Cpg*V0pg + Cmet*V0met - 1 (2)

% 0 = VfCpo*Cpo - Vt*Cpo - k*Cpo*V (3)
% 0 = VfCh20*Ch2o - Vt*Ch2o - k*Cpo*V (4)
% 0 = - Vt*Cpg + k*Cpo*V (5)
% 0 = VfCmet*Cmet - Vt*Cmet  (6)

result = fsolve(@(x)prob1(x,mass,MolWtPO,MolWtH2O,MolWtPG,MolWtMet,molVolPO,molVolH2O,molVolPG,molVolMet,molFeedPO,molFeedH2O,molFeedMet,k),[1 1 1 1 1000 7000]);

Cpo = result(1); %0.671236063161155 Mol/L
Ch2o = result(2); % 45.570908455934610 Mol/L
Cpg = result(3); %1.816375561125502 Mol/L
Cmet = result(4); % 8.999628129348908e-24 Mol/L
Vol = result(5); % 792.6543716737496 L
VolFlow = result(6); % 7729.071607055675 L/h

Cpo_0 = (-VolFlow*Cpo-k*Cpo*Vol)/VfPO;

%Conversion of PO
XPO = (Cpo_0-Cpo)/Cpo_0; %104.5384622291502 %
%production rate
RCPO = k*Cpo; %17.711246261848963 mol / L * h

%it seams that the best choice would be constant volume, as its production
%rate is quite larger than at constant mass, our volumetric flow is larger
%but it is less efficient

%(b)

%assuming all densities are equal to water
densH2O = 1.000;   % g/cm^3
densPO = densH2O;  % g/cm^3
densPG = densH2O;  % g/cm^3
densMet = densH2O; % g/cm^3

%molar volume of the components
molVolPO = (MolWtPO / densPO)/1000; % L / mol
molVolH2O = (MolWtH2O / densH2O)/1000; % L / mol
molVolPG = (MolWtPG / densPG)/1000; % L / mol
molVolMet = (MolWtMet/densMet)/1000; % L / mol

%molar feed rate of each component
molFeedPO = VfPO / molVolPO; %Mol / h
molFeedH2O = VfH2O/molVolH2O; %Mol / h
molFeedPG = VfPG/molVolPG ; %Mol / h
molFeedMet =  VfMet/molVolMet; %Mol / h

%inital concentration of each component
cPO_0 = molFeedPO / VfPO; %mol / L
cH2O_0 = molFeedH2O / VfH2O; %mol / L
cPG_0 = 0;%molFeedPG / VfPG;(NaN) %mol / L
cMet_0 = molFeedMet / VfMet; %mol / L

%concentration of each component at the reactor exit (Ca)
cPO = cPO_0 / (1 + k * (V/VfT)); % 4.451248457264060 mol / L
cH2O = cH2O_0 / (1 + k * (V/VfT)); % 14.346754184123006 mol / L
cPG = cPO_0 - cPO; %12.766382396730428 mol /L

%Conversion of PO
xPO = (cPO_0-cPO)/cPO_0; % 0.741471489602103
%Production rate
rCPO = k*cPO; %117.4507180499200 mol / L * h

%larger production rate when compared to using the "correct values"
%with a error of 17% in the production rate, would be best to use the more
%accurate values

%Problem 2
%plot concenration vs time of each component
%find time to max B

% A -> B <-> C

% Species  | Feed |        Exit (Ni)     |          yi          |     
%---------------------------------------------------------------|
% A        | Ca0  |        Ca0 - z1      |  (Ca0 - z1)/(Ca0)    |
% B        |  0   |    0 + z1 - z2 + z3  | (z1 - z2 + z3)/(Ca0) |
% C        |  0   |       0 + z2 - z3    |   (z2 - z3)/(Ca0)    |
%----------|------|----------------------|----------------------|
% Total    | Ca0  |         Ca0          |          1           |

Ca0 = 1.0;
Cb0 = 0.0;
Cc0 = 0.0;
k1 = 0.5; %min-1
k2 = 0.25; %min-1
k3 = 0.1; %min-1

t_range = linspace(0,20,120);
[y,x] = ode45(@(t,x)prob2(t,x,k1,k2,k3),t_range,[Cb0 Cc0 Ca0]);
Cb = x(:,1);
Cc = x(:,2);
Ca = x(:,3);

maxCb = max(Cb);
idxMaxCb = find(Cb==maxCb);

hold on
plot(t_range,Ca,"- k") %ploting of Ca
plot(t_range,Cb,"-- k") %ploting of Cb
plot(t_range,Cc,": k") %ploting of Cc
plot(t_range(idxMaxCb),maxCb,"d k","MarkerFaceColor","k"); %ploting of the max Cb
xlabel("time (min)")
ylabel("Concetration - Cj")
title("Concentration of components over time")
legend("Ca","Cb","Cc","Max Cb = 0.5242")
hold off

conc = fsolve(@(x)prob2_1(x,k1,k2,k3),[0.2 0.7 0.1]); %eq prediction
CaPred = conc(3); %-1.659348134926647e-06
CbPred = conc(1); %0.268961096216199
CcPred = conc(2); %0.672415561462095

%as it can be noted the Eq prediction does aproximate the value of the
%concentration of the components when time increases, in t=20 minutes the
%diference was minimal, and it could be seen as a a accurate decription of
%the reaction at equilibrium.

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

alfaCO2 = (1 + kCO2*(1-sqrt(T/TcCO2))).^2;
alfaCH4 = (1 + kCH4*(1-sqrt(T/TcCH4))).^2;

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

fulgacityCO2 = phiCO2 * yCO2 * Pt; %f = 208.7810648391098 bar
fulgacityCH4 = phiCH4 * yCH4 * Pt; %f = 265.5299366509511 bar

function f = cmp(z,A,B)
f = z.^3 - (1-B)*z.^2 + (A-3*(B.^2)-2*B)*z - (A*B -B.^2 - B.^3);
end
function w = prob1(x,Mass,Mpo,Mh2o,Mpg,Mmet,V0po,V0h2o,V0pg,V0met,MolFeedPO,MolFeedH2O,MolFeedMet,k)

w(1)= (Mass/(x(1)*Mpo+x(2)*Mh2o+x(3)*Mpg+x(4)*Mmet)) - x(5);
w(2)= x(1)*V0po + x(2)*V0h2o + x(3)*V0pg + x(4)*V0met - 1;
w(3) = MolFeedPO - x(6)*x(1) - k*x(1)*x(5);
w(4) = MolFeedH2O - x(6)*x(2) - k*x(1)*x(5);
w(5) = - x(6)*x(3) + k*x(1)*x(5);
w(6) = MolFeedMet*x(4) - x(6)*x(4);

end
function u = prob2(t,x,k1,k2,k3)
    dcadt = -k1*x(3);
    dcbdt = k1*x(3) - k2*x(1) + k3*x(2);
    dccdt = k2*x(1) - k3*x(2);
    u = [dcbdt;dccdt;dcadt];
end
function u = prob2_1(x,k1,k2,k3)
    u(1) = -k1*x(3);
    u(2) = k1*x(3) - k2*x(1) + k3*x(2);
    u(3) = k2*x(1) - k3*x(2);
end