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

%stoic table
% Species | Fi0  | Exit (Fi) |     Exit (Ci)     |     
%---------|------|-----------|-------------------|
% A       |  F0  |    F0 - z |   F0-F0*Xa / v    |
% B       |   0  |    0.3z   |   0.3*F0*Xa / v   |
% C       |   0  |    0.065z |  0.065*F0*Xa / v  |
% D       |   0  |    0.6675z| 0.6675*F0*Xa / v  |
% E       |   0  |    0.653z |  0.635*F0*Xa / v  |
% F       |   0  |    0.3z   |   0.3*F0*Xa / v   |
%---------|------|-----------|-------------------|
% Total   |  F0  |    Ft     |  F0-0.9855*F0*Xa/v|

% Fa0 = F0
% Fa = F0 - F0*Xa
% Fa = F0 - z
% then
% z = F0*Xa (extent of reaction)

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
Q = Q/3600; %cal/(s*m2)

%pipes data
%2 inch diameter pipe
di_2 = 2.07; %in        Internal diameter
ai_2 = 3.36; %in2       Internal Transverse Area
do_2 = 2.375; %in       External Diameter
wpft_2 = 3.65; %lb/ft   Weight per length
di_2 = di_2 / 39.37;%m
do_2 = do_2 / 39.37;%m
ai_2 = ai_2 / 1550;%m2
ff_2 = 0.0050;%         Friction factor
%4 inch diameter pipe
di_4 = 4.03; %in        Internal diameter
ai_4 = 12.73;%in2       Internal Transverse Area
do_4 = 4.500; %in       External Diameter
wpft_4 = 10.79;%lb/ft   Weight per length
di_4 = di_4 / 39.37;%m
do_4 = do_4 / 39.37;%m
ai_4 = ai_4 / 1550;%m2
ff_4 = 0.0044;%         Friction factor
%6 inch diameter pipe
di_6 = 6.07; %in        Internal diameter
ai_6 = 28.89;%in2       Internal Transverse Area
do_6 = 6.625; %in       External Diameter
wpft_6 = 18.97;%lb/ft   Weight per length
di_6 = di_6 / 39.37;%m
do_6 = do_6 / 39.37;%m
ai_6 = ai_6 / 1550;%m2
ff_6 = 0.0044;%         Friction factor

%heats of reactions
%gaussian calculated (298 K)
dH_H2_gau = -1.162033; %Hartrees
dH_CH4_gau = -40.469367; %hartrees
dH_C3H8_gau = -119.034669; %hartrees
dH_C3H6_gau = -117.822454; %hartrees
dH_C2H6_gau = -79.750776; %hartress
dH_C2H4_gau = -78.552125; %hartress

dH_prod = vB*dH_C3H6_gau + vC*dH_C2H6_gau + vD*dH_C2H4_gau + vE*dH_CH4_gau + vF*dH_H2_gau; %hartrees
dH_reac = vA*dH_C3H8_gau; %hartrees

dH_Gau_298 = dH_prod - dH_reac; %hartrees
dH_Gau_298 = dH_Gau_298 * 627.5095; %Kcal/mol
dH_Gau_298 = dH_Gau_298 * 1000; %cal/mols

%length of the reactor (m)
len_span = linspace(0,426,400); %426 m ~ 1400 ft

%initial conditions
y0 = [F0,0,0,0,0,0,T0,P0];

%RUNS%
%Gaussian dH
%2 in
Tr = 298; %K
dH_ref = dH_Gau_298; %cal/mol
%Pipe data:
A = ai_2;
D = di_2;
DO = do_2;
ff = ff_2;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_2_in = x(:,7); %K
P_2_in = x(:,8); %Pa
P_2_in = P_2_in./6895;%psia
Fa_2_in = x(:,1);%mol/s
Xa_2_in = (F0-Fa_2_in)./F0;
Xa_2_in = Xa_2_in*100; % (%)
L_2_in = l * 3.281; % ft
%Plots
figure;
subplot(2,2,1);
plot(L_2_in,T_2_in);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_2_in,P_2_in);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_2_in,Xa_2_in);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("2 inches diameter pipe");

%4 in
Tr = 298; %K
dH_ref = dH_Gau_298; %cal/mol
%Pipe data:
A = ai_4;
D = di_4;
DO = do_4;
ff = ff_4;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_4_in = x(:,7); %K
P_4_in = x(:,8); %Pa
P_4_in = P_4_in./6895;%psia
Fa_4_in = x(:,1);%mol/s
Xa_4_in = (F0-Fa_4_in)./F0;
Xa_4_in = Xa_4_in*100; % (%)
L_4_in = l * 3.281; % ft
%Plots
figure;
subplot(2,2,1);
plot(L_4_in,T_4_in);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_4_in,P_4_in);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_4_in,Xa_4_in);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("4 inches diameter pipe");

%6 in
Tr = 298; %K
dH_ref = dH_Gau_298; %cal/mol
%Pipe data:
A = ai_6;
D = di_6;
DO = do_6;
ff = ff_6;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_6_in = x(:,7); %K
P_6_in = x(:,8); %Pa
P_6_in = P_6_in./6895;%psia
Fa_6_in = x(:,1);%mol/s
Xa_6_in = (F0-Fa_6_in)./F0;
Xa_6_in = Xa_6_in*100; % (%)
L_6_in = l * 3.281; %ft
%Plots
figure;
subplot(2,2,1);
plot(L_6_in,T_6_in);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_6_in,P_6_in);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_6_in,Xa_6_in);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("6 inches diameter pipe");

%grouped plots
figure;
subplot(2,2,1);
hold on;
plot(L_2_in,T_2_in,"k-");
plot(L_4_in,T_4_in,"k--");
plot(L_6_in,T_6_in,"k:");
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
legend("2 in","4 in","6 in");
hold off;
subplot(2,2,2);
hold on;
plot(L_2_in,P_2_in,"k-");
plot(L_4_in,P_4_in,"k--");
plot(L_6_in,P_6_in,"k:");
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
legend("2 in","4 in","6 in");
hold off;
subplot(2,2,[3,4]);
hold on;
plot(L_2_in,Xa_2_in,"k-");
plot(L_4_in,Xa_4_in,"k--");
plot(L_6_in,Xa_6_in,"k:");
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
legend("2 in","4 in","6 in");
hold off;
sgtitle("Grouped plots");

%6 in  - DH exam
Tr = 8.664833333333333e2; %K
dH_ref = 21.96*1000; %cal/mol
%Pipe data:
A = ai_6;
D = di_6;
DO = do_6;
ff = ff_6;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_6_in_dh = x(:,7); %K
P_6_in_dh = x(:,8); %Pa
P_6_in_dh = P_6_in_dh./6895;%psia
Fa_6_in_dh = x(:,1);%mol/s
Xa_6_in_dh = (F0-Fa_6_in_dh)./F0;
Xa_6_in_dh = Xa_6_in_dh*100; % (%)
L_6_in_dh = l * 3.281; %ft
%Plots
figure;
subplot(2,2,1);
plot(L_6_in_dh,T_6_in_dh);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_6_in_dh,P_6_in_dh);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_6_in_dh,Xa_6_in_dh);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("6 inches diameter pipe (with exan given dH)");

%comparison between the gaussian and the given value
figure;
subplot(3,1,1);
hold on;
plot(L_6_in,Xa_6_in,"k-");
plot(L_6_in_dh,Xa_6_in_dh,"k--");
title("Conversion vs Length (6 in)");
xlabel("Length (ft)");
ylabel("Conversion (%)");
legend("Gaussian","Given");
hold off;
subplot(3,1,2);
hold on;
plot(L_4_in,Xa_4_in,"k-");
plot(L_4_in_dh,Xa_4_in_dh,"k--");
title("Conversion vs Length (4 in)");
xlabel("Length (ft)");
ylabel("Conversion (%)");
legend("Gaussian","Given");
hold off;
subplot(3,1,3);
hold on;
plot(L_2_in,Xa_2_in,"k-");
plot(L_2_in_dh,Xa_2_in_dh,"k--");
title("Conversion vs Length (2 in)");
xlabel("Length (ft)");
ylabel("Conversion (%)");
legend("Gaussian","Given");
hold off;
sgtitle("Comparison of given Heat of reaction vs Gaussian calculated");



%4 in  - DH exam
Tr = 8.664833333333333e2; %K
dH_ref = 21.96*1000; %cal/mol
%Pipe data:
A = ai_4;
D = di_4;
DO = do_4;
ff = ff_4;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_4_in_dh = x(:,7); %K
P_4_in_dh = x(:,8); %Pa
P_4_in_dh = P_4_in_dh./6895;%psia
Fa_4_in_dh = x(:,1);%mol/s
Xa_4_in_dh = (F0-Fa_4_in_dh)./F0;
Xa_4_in_dh = Xa_4_in_dh*100; % (%)
L_4_in_dh = l * 3.281; % ft
%Plots
figure;
subplot(2,2,1);
plot(L_4_in_dh,T_4_in_dh);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_4_in_dh,P_4_in_dh);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_4_in_dh,Xa_4_in_dh);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("4 inches diameter pipe (with exam given dH)");

%2 in  - DH exam
Tr = 8.664833333333333e2; %K
dH_ref = 21.96*1000; %cal/mol
%Pipe data:
A = ai_2;
D = di_2;
DO = do_2;
ff = ff_2;
%ODE run
[l,x] = ode15s(@(L,x)PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO),len_span,y0);
T_2_in_dh = x(:,7); %K
P_2_in_dh = x(:,8); %Pa
P_2_in_dh = P_2_in_dh./6895;%psia
Fa_2_in_dh = x(:,1);%mol/s
Xa_2_in_dh = (F0-Fa_2_in_dh)./F0;
Xa_2_in_dh = Xa_2_in_dh*100; % (%)
L_2_in_dh = l * 3.281; % ft
%Plots
figure;
subplot(2,2,1);
plot(L_2_in_dh,T_2_in_dh);
title("Temperature vs Length");
xlabel("Length (ft)");
ylabel("Temperature (K)");
subplot(2,2,2);
plot(L_2_in_dh,P_2_in_dh);
title("Pressure vs Length");
xlabel("Length (ft)");
ylabel("Pressure (psia)");
subplot(2,2,[3,4]);
plot(L_2_in_dh,Xa_2_in_dh);
title("Conversion vs Length");
xlabel("Length (ft)");
ylabel("Conversion (%)");
sgtitle("2 inches diameter pipe (with exam given dH)");

%PFR Function
function f=PFR(L,x,T0,P0,F0,M0,v0,Q,Tr,dH_ref,A,ff,D,DO)
%variables
Fa = x(1);
Fb = x(2);
Fc = x(3);
Fd = x(4);
Fe = x(5);
Ff = x(6);
T = x(7);
P = x(8);
%R gas constant
R = 1.98720425864083;%cal/(mol*K)
%rate constant (s-1)
k = 3.98e12*exp(-59100/(R*T));
%conversion
Xa = (F0-Fa)/F0;
%total molar flow
Ft = Fa+Fb+Fc+Fd+Fe+Ff; %mol/s
%volumetric flow
v = v0*(Ft/F0)*(T/T0)*(P0/P); %m3/s
%concentrations
Ca = (F0-F0*Xa)/v; %mol/m3
Cb = (0.3*F0*Xa) / v; %mol/m3
Cc = (0.065*F0*Xa) / v; %mol/m3
Cd = (0.6675*F0*Xa) / v; %mol/m3
Ce = (0.635*F0*Xa) / v; %mol/m3
Cf = (0.3*F0*Xa) / v; %mol/m3
%reaction rate
ra = -k*Ca; %mol/s*m3
rb = k*Cb; %mol/s*m3
rc = k*Cc; %mol/s*m3
rd = k*Cd; %mol/s*m3
re = k*Ce; %mol/s*m3
rf = k*Cf; %mol/s*m3
%Density 
rho = M0/v; %Kg/m3
%heat of reaction
dH = dH_ref + integral(@dCP,Tr,T);
%Heat capacities (cal/molK)
cpA = 21.14 + 0.02056*T;
cpB = 17.88 + 0.01645*T;
cpC = 13.34 + 0.01589*T;
cpD = 12.29 + 0.01022*T;
cpE = 6.98 + 0.01012*T;
cpF = 6.42 + 0.00082*T;
%Temperature along the length of the reactor
FiCpi = Fa*cpA + Fb*cpB + Fc*cpC + Fd*cpD + Fe*cpE + Fd*cpD + Ff*cpF;
dTdl = (ra*dH - Q*(1/(2*pi*(DO/2))))/(FiCpi) * A;
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
G = M0 / A;
dPdl = -(2*ff*(G^2))/(rho*D);
%final
f = [dFadl;dFbdl;dFcdl;dFddl;dFedl;dFfdl;dTdl;dPdl];
end

function x=dCP(T)
    %heat Capacities (cal/mol*K)
    cpA = 21.14 + 0.02056*T;
    cpB = 17.88 + 0.01645*T;
    cpC = 13.34 + 0.01589*T;
    cpD = 12.29 + 0.01022*T;
    cpE = 6.98 + 0.01012*T;
    cpF = 6.42 + 0.00082*T;
    dCp = -cpA + 0.3*cpB + 0.065*cpC + 0.6675*cpD + 0.653*cpE + 0.3*cpF;
    x = dCp;
end



