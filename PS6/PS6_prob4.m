%Problem set 6
%problem 4

%RAW 6.17

% A - Cl2
% B - C3H6
% C - C3H5Cl
% D - HCl
% E - 12C3H6Cl2

% A + B -> C + D  (1)
% A + B -> E      (2)

% Species | Fi0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       | Fa0  |    Fa     |Fa0-z1-z2 |
% B       | Fb0  |    Fb     |Fb0-z1-z2 |
% C       |   0  |    Fc     |    z1    |
% D       |   0  |    Fd     |    z1    |
% E       |   0  |    Fe     |    z2    |
%---------|------|-----------|----------|
% Total   | Ft0  |    Ft     | Ft0 - z2 |


HfA = 0.0; %Kcal/Mol
HfB = 4.88; %Kcal/Mol
HfC = -0.15; %Kcal/Mol
HfD = -22.06; %Kcal/Mol
HfE = -39.60; %Kcal/Mol

aA = 6.432;
aB = 0.866;
aC = 0.604;
aD = 7.235;
aE = 2.496;

bA = 0.8082 / 1e2;
bB = 5.602 / 1e2;
bC = 7.277 / 1e2;
bD = -0.172 / 1e2;
bE = 8.729 / 1e2;

cA = -0.9241 / 1e5;
cB = -2.771 / 1e5;
cC = -5.442 / 1e5;
cD = 0.2976 / 1e5;
cE = -6.219 / 1e5;

dA = 3.695 / 1e9;
dB = 5.266 / 1e9;
dC = 17.42 / 1e9;
dD = -0.931 / 1e9;
dE = 18.49 / 1e9;

dHr1_0 = -HfA - HfB + HfC + HfD; %Kcal/Mol
dHr2_0 = -HfA - HfB + HfE; %Kcal/Mol

dHr1_0 = dHr1_0 * 9.223e18; %Kcal/mol to BTU/lbmol
dHr2_0 = dHr2_0 * 9.223e18; %Kcal/mol to BTU/lbmol

Pt0 = 2; %atm
len = 25; %ft
diameter = 2; %in
T0 = 392; %F
Ft0 = 0.85;% Lbmol/hr

diameter = diameter / 12; %in to ft
Ac = (pi * diameter^2) / 4; %ft2

ap = diameter*pi;

Ta = 386.6 + 459.67; %R

Vr = Ac * len; %ft3
T0 = T0 + 459.67; % R
R = 0.730240507295273; %atm⋅ft3⋅lbmol−1⋅°R−1
V0 = (Ft0 * R * T0)/Pt0 ; %ft3/hr

yA0 = 1/5;
yB0 = 4/5;

Fa0 = Ft0 * yA0; %lbmol/hr
Fb0 = Ft0 * yB0; %lbmol/hr

U = 5; %Btu/hr Ft2 F


L_span = linspace(0,len,500);
y0 = [Fa0 Fb0 0.0 0.0 0.0 T0];
[L,x] = ode15s(@(L,x)PFR(L,x,Ac,V0,T0,R,dHr1_0,U,Ta,ap,dHr2_0,aA,aB,aC,aD,aE,bA,bB,bC,bD,bE,cA,cB,cC,cD,cE,dA,dB,dC,dD,dE),L_span,y0);

Fa = x(:,1);
Fb = x(:,2);
Fc = x(:,3);
Fd = x(:,4);
Fe = x(:,5);
T = x(:,6);

hold on
plot(L,Fa);
plot(L,Fc);
plot(L,Fe);
legend("Cl_2","C_3H_5Cl","C_3H_6Cl_2");
xlabel("Length (ft)");
ylabel("lbmolflow (lbmol/h)");
title("Molar Flow versus length");
hold off

function f = PFR(L,x,Ac,V0,T0,R,dHr1_0,U,Ta,ap,dHr2_0,aA,aB,aC,aD,aE,bA,bB,bC,bD,bE,cA,cB,cC,cD,cE,dA,dB,dC,dD,dE)

FA = x(1);
FB = x(2);
FC = x(3);
FD = x(4);
FE = x(5);
T = x(6);

k1 = 2.06e5 * exp(-27200/(R*T));
k2 = 11.7 * exp(-6860/(R*T));

V = V0 * (T/T0);

Pa = (FA * R * T)/V;
Pb = (FB * R * T)/V;

r1 = k1 * Pa * Pb;
r2 = k2 * Pa * Pb;

ra = -r1 -r2;
rb = -r1 -r2;
rc = r1;
rd = r1;
re = r2;

Cpa = (aA + bA*T + cA*T^2 + dA*T^3);
Cpb = (aB + bB*T + cB*T^2 + dB*T^3);
Cpc = (aC + bC*T + cC*T^2 + dC*T^3);
Cpd = (aD + bD*T + cD*T^2 + dD*T^3);
Cpe = (aE + bE*T + cE*T^2 + dE*T^3);

dCp1 = @(T) (-Cpa-Cpb+Cpc+Cpd);
dCp2 = @(T) (-Cpa-Cpb+Cpe);

deltaH1 = (integral(dCp1,T0,T,'ArrayValued',true))/1000; %Kcal/mol
deltaH2 = (integral(dCp2,T0,T,'ArrayValued',true))/1000; %Kcal/mol

dHr1 = dHr1_0 + (deltaH1 * 9.223e18);
dHr2 = dHr2_0 + (deltaH2 * 9.223e18);

%adiabatic Case
Q = 0;

%non adiabatic
% Q = U*ap*(Ta-T);

dFAdL = ra * Ac;
dFBdL = rb * Ac;
dFCdL = rc * Ac;
dFDdL = rd * Ac;
dFEdL = re * Ac;
dTdL = (Q + ((-dHr1)*r1+(-dHr2)*r2)*Ac)/ (FA * Cpa + FB * Cpb);
f = [dFAdL;dFBdL;dFCdL;dFDdL;dFEdL;dTdL];

end
