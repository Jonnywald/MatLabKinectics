%Final Exam%
%Guilherme Bertola%

% C3H8 -> 0.3 C3H6 + 0.065 C2H6 + 0.6675 C2H4 + 0.635 CH4 + 0.3 H2
% A -> 0.3 B + 0.065 C + 0.6675 D + 0.635 E + 0.3 F

%Letter | Molecule
%  A    |   C3H8
%  B    |   C3H6
%  C    |   C2H6
%  D    |   C2H4
%  E    |   CH4
%  F    |   H2

%stoic indexes
vA = 1;
vB = 0.3;
vC = 0.065;
vD = 0.6675;
vE = 0.635;
vF = 0.3;

%known parameters
%conversion
Xa = 0.8;
%Initial Temperature
T0 = 1100; %F
T0 = (T0-32) * 5/9 + 273.15; %K
%Inital Pressure
P0 = 60;% psia
P0 = P0 * 6894.76; %Pa
%molar mass of propane
Ma = 44.097;%g/mol
Ma = Ma/1000; %kg/mol
%Feed Rate (mass)
M0 = 7000;% lb/h
M0 = M0/2.205; %Kg/h
M0 = M0/3600; %Kg/s
%Feed Rate (mol)
F0 = M0/Ma; %mol/s
%R gas constant
R = 8.31446261815324; %m3⋅Pa⋅K−1⋅mol−1
%concentration
c = P0/(R*T0); %mol/m3
%volumetric flow
v0 = F0/c; %m3/s
%Heat flux
Q = 2.52e6; %cal/(h*ft2)
Q = Q*10.764; %cal/(h*m2)

%pipes data
%2 inch diameter pipe
di_2 = 2.07; %in        Internal diameter
ai_2 = 3.36; %in2       Internal Transverse Area
di_2 = di_2 / 39.37;%m
ai_2 = ai_2 / 1550;%m2
%4 inch diameter pipe
di_4 = 4.03; %in        Internal diameter
ai_4 = 12.73;%in2       Internal Transverse Area
di_4 = di_4 / 39.37;%m
ai_4 = ai_4 / 1550;%m2
%6 inch diameter pipe
di_6 = 6.07; %in        Internal diameter
ai_6 = 28.89;%in2       Internal Transverse Area
di_6 = di_6 / 39.37;%m
ai_6 = ai_6 / 1550;%m2
%heats of reactions
%gaussian calculated (298 K)
dH_H2_gau = -1.162033; %Hartrees
dH_CH4_gau = -40.469367; %hartrees
dH_C3H8_gau = -119.034669; %hartrees
dH_C3H6_gau = -117.820119; %hartrees ** redo not converged
dH_C2H6_gau = -79.750776; %hartress
dH_C2H4_gau = -78.552125; %hartress

dH_prod = vB*dH_C3H6_gau + vC*dH_C2H6_gau + vD*dH_C2H4_gau + vE*dH_CH4_gau + vF*dH_H2_gau; %hartrees
dH_reac = vA*dH_C3H8_gau; %hartrees

dH_Gau_298 = dH_prod - dH_reac; %hartrees
dH_Gau_298 = dH_Gau_298 * 627.5095; %Kcal/mol


function f=PFR(L,x,F0,M0,v0,Q,Tr,dH_ref,A,ff,D)
%variables
Fa = x(1);
Fb = x(2);
Fc = x(3);
Fd = x(4);
Fe = x(5);
Ff = x(6);
T = x(7);
P = x(8);
%Density of Propane gas
rho = 1.898; %Kg/m3
%R gas constant
R = 1.98720425864083;%cal/(mol*K)
%rate constant (s-1)
k = 3.98e12*exp(-59100/(R*T));
%heat Capacities (cal/mol*K)
cpA = 21.14 + 0.02056*T;
cpB = 17.88 + 0.01645*T;
cpC = 13.34 + 0.01589*T;
cpD = 12.29 + 0.01022*T;
cpE = 6.98 + 0.01012*T;
cpF = 6.42 + 0.00082*T;
dCp = -cpA + 0.3*cpB + 0.065*cpC + 0.6675*cpD + 0.653*cpE + 0.3*cpF; 
%total molar flow
Ft = Fa+Fb+Fc+Fd+Fe+Ff; %mol/s
%volumetric flow
v = v0*(Ft/F0)*(T/T0)*(P0/P); %m3/s
%concentrations
Ca = Fa/v; %mol/m3
Ca0 = F0/v0; %mol/m3
Cb = (Ca0 - Ca) * 0.3; %mol/m3
Cc = (Ca0 - Ca) * 0.065; %mol/m3
Cd = (Ca0 - Ca) * 0.6675; %mol/m3
Ce = (Ca0 - Ca) * 0.653; %mol/m3
Cf = (Ca0 - Ca) * 0.3; %mol/m3
%reaction rate
ra = -k*Ca; %mol/s*m3
rb = k*Cb; %mol/s*m3
rc = k*Cc; %mol/s*m3
rd = k*Cd; %mol/s*m3
re = k*Ce; %mol/s*m3
rf = k*Cf; %mol/s*m3
%heat of reaction
dH = dH_ref + dCp*(T-Tr);
%Temperature along the length of the reactor
FiCpi = Fa*cpA + Fb*cpB + Fc*cpC + Fd*cpD + Fe*cpE + Fd*cpD + Ff*cpF;
dTdl = (ra*dH - Q)/(FiCpi) * A;
%Molar flow of A along the length of the reactor
dFadl = ra * A;
%Molar flow of B along the length of the reactor
dFbdl = rb * A;
%Molar flow of C along the length of the reactor
dFcdl = rc * A;
%Molar flow of D along the length of the reactor
dFddl = rd * A;
%Molar flow of E along the length of the reactor
dFedl = re * A;
%Molar flow of F along the length of the reactor
dFfdl = rf * A;
%pressure drop along the length of the reactor
G = M0 * A;
dPdl = -(2*ff*G^2)/(rho*D);
%final
f = [dFadl;dFbdl;dFcdl;dFddl;dFedl;dFfdl;dTdl;dPdl];
end



