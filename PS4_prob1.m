%problem set 4
%Problem 1

%ethane pyrolysis

% C2H6 + NO <-> C2H5 + HNO  (1)
% C2H5 -> H + C2H4          (2)
% H + C2H6 -> C2H5 + H2     (3)
% H + NO <->  HNO           (4)

% Species | MolFeed Fi0 | Exit (Fi) |     yi    |     
%---------|-------------|-----------|-----------|
% C2H6    |   FC2H60    |   FC2H6   |  FC2H6/FT |
% NO      |    FNO0     |   FNO     |  FNO/FT   |
% C2H5    |     0       |   FC2H5   |  FC2H5/FT |
% HNO     |     0       |   FHNO    |  FHNO/FT  |
% H       |     0       |   FH      |  FH/FT    |
% C2H4    |     0       |   FC2H4   |  FC2H4/FT |
% H2      |     0       |   FH2     |  FH2/FT   |
%---------|-------------|-----------|-----------|
% Total   |FC2H60 + FNO0|    FT     |      1    |

% r1 = k1*cC2H6*cNO - kr1*cC2H5*cHNO
% r2 = k2*cC2H5
% r3 = k3*cH*cC2H6
% r4 = k4*cH*cNO - kr4*cHNO

A1 = 1.0e14;    %cm3/mol*s
Ar1 = 1.0e12;   %cm3/mol*s
A2 = 3.0e14;    %s-1
A3 = 3.4e12;    %cm3/mol*s
A4 = 1.0e12;    %cm3/mol*s
Ar4 = 1.0e13;   %s-1

E1 = 217.6;     %KJ/mol
Er1 = 0.0;      %KJ/mol
E2 = 165.3;     %KJ/mol
E3 = 28.5;      %KJ/mol
E4 = 0.0;       %KJ/mol
Er4 = 200.8;    %KJ/mol

V = 1500;       %cm3
P = 1.0;        %atm
T = 1050;       %K
V0 = 600;       %cm3/s
R_KJ = 8.31446261815324e-3; % KJ / K*Mol
R_atm = 82.0573660809596; % cm3⋅atm⋅K−1⋅mol−1

k1 = A1*exp(-E1/(R_KJ*T));
kr1 = Ar1*exp(-Er1/(R_KJ*T));
k2 = A2*exp(-E2/(R_KJ*T));
k3 = A3*exp(-E3/(R_KJ*T));
k4 = A4*exp(-E4/(R_KJ*T));
kr4 = Ar4*exp(-Er4/(R_KJ*T));

yNO = 5e-6;     %5 ppm
yC2H6 = 1 - yNO;

FNO0 = yNO*(P/(R_atm*T))*V0; %3.481884678486923e-08 mol/s
FC2H60 = yC2H6*(P/(R_atm*T))*V0; %0.006963734538127 mol/s

cC2H60 = FC2H60/V0; %1.160622423021177e-05 mol/cm3
cNO0 = FNO0/V0; %5.803141130811539e-11 mol/cm3

v_span = linspace(0,V,500);

[v,c] = ode15s(@(v,c)func(v,c,k1,kr1,k2,k3,k4,kr4,V0),v_span,[cC2H60 cNO0 0 0 0 0 0]);

cC2H6 = c(:,1);
cNO = c(:,2);
cC2H5 = c(:,3);
cHNO = c(:,4);
cH = c(:,5);
cC2H4 = c(:,6);
cH2 = c(:,7);

%conversion calculation
XC2H6 = (cC2H60-cC2H6)./cC2H60;

%doing the same with 5% NO
yNO = 0.05;
yC2H6 = 1 - yNO;
FNO0 = yNO*(P/(R_atm*T))*V0;
FC2H60 = yC2H6*(P/(R_atm*T))*V0;
cC2H60 = FC2H60/V0;
cNO0 = FNO0/V0;
[v,c] = ode15s(@(v,c)func(v,c,k1,kr1,k2,k3,k4,kr4,V0),v_span,[cC2H60 cNO0 0 0 0 0 0]);
cC2H6_p = c(:,1);
XC2H6_p = (cC2H60-cC2H6_p)./cC2H60;

hold on
plot(v_span,XC2H6_p,"-");
plot(v_span,XC2H6,"--");
xlabel("Volume (cm3)");
ylabel("Conversion");
title("Conversion by volume of PFR");
legend("XC_2H_6 - 5% NO","XC_2H_6 - 5ppm NO");
hold off

%reducing the concentration of NO in the feed to 5 ppm, lowers the
%conversion oh ethane siginificantly (from what as about 80% to 35%) as expected,
%because NO is a Key component in the ethane pyrolysis

function f = func(v,c,k1,kr1,k2,k3,k4,kr4,V0)

r1 = k1*c(1)*c(2) - kr1*c(3)*c(4);
r2 = k2*c(3);
r3 = k3*c(5)*c(1);
r4 = k4*c(5)*c(2) - kr4*c(4);

dcC2H6dv = (-r1 -r3)/V0;
dcNOdv = (-r1 -r4)/V0;
dcC2H5dv = (r1 -r2 + r3)/V0;
dcHNOdv = (r1 + r4)/V0;
dcHdv = (r2 - r3 - r4)/V0;
dcC2H4dv = r2/V0;
dcH2dv = r3/V0;

f = [dcC2H6dv;dcNOdv;dcC2H5dv;dcHNOdv;dcHdv;dcC2H4dv;dcH2dv];

end