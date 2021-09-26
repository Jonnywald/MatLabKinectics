%problem set 4
%problem 4
% ammonia Oxidation

% Species | MolFeed Fi0| Exit (Fi) |    yi    |     
%---------|------------|-----------|----------|
% NH3     |    FNH30   |   FNH3    |  FNH3/FT |
% O2      |    FO20    |    FO2    |  FO2/FT  |
% NO      |     0      |    FNO    |  FNO/FT  |
% H2O     |     0      |   FH2O    |  FH2O/FT |
% N2      |     0      |    FN2    |  FN2/FT  |
% NO2     |     0      |   FNO2    |  FNO2/FT |
%---------|------------|-----------|----------|
% Total   |FNH30 + FO20|    FT     |     1    |

StoichiMat = [-4 -5 4 6 0 0;
              -2 -1.5 0 3 1 0;
              0 -1 -2 0 0 2;
              -4 0 -6 6 5 0]; %stoichiometric matrix

% -R1Nh3 = k1*Cnh3*Co2^2
% -R2Nh3 = k2*Cnh3*Co2
% -R3O2 = k3*Cno2^2*Co2
% -R4NO = k4*Cn0*Cnh3^2/3

V0 = 10; %Feed rate (L/min)
V = 10; %volume of the reactor (L)

k1 = 5.0; %(M^3/kmol)^2 / min
k2 = 2.0; %(M^3/kmol-min)
k3 = 10.0;%(M^3/kmol)^2 / min
k4 = 5.0; %(M^3/kmol)^2/3 / min
% the rates need to be converted to the same units used int the problem
% M3 -> L , Kmol -> mol

k1 = k1 * 1000000; % conversion between M3 to Liters (squared)
k2 = k2 * 1000;    % conversion between M3 to Liters
k3 = k3 * 1000000; % conversion between M3 to Liters (squared)
k4 = k4 * 100;     % conversion between M3 to Liters (to the power of 2/3)

k1 = k1 * (1/1000000); % conversion between Kmol to mol (squared)
k2 = k2 * (1/1000);    % conversion between Kmol to mol
k3 = k3 * (1/1000000); % conversion between Kmol to mol (squared)
k4 = k4 * (1/100);     % conversion between Kmol to mol (to the power of 2/3)

% the values stayed the same in the end


cNH30 = 1.0; % Initial Concentration of NH3 (mol/L)
cO20 = 1.0; % Initial Concentration of O2 (mol/L)
FNH30 = cNH30 * V0; %Initial Mol Flow of NH3 = 10 mol/min
FO2 = cO20 * V0; %Initial Mol Flow of O2 = 10 mol/min

% dcNh3dv = rnh3 = -r1 -r2 -r4 / V0
% dcO2dv = rO2 = -r1 -r2 -r3 / V0
% dcnodv = rno = r1 -r3 -r4 / V0
% dcH2Odv = rh2o = r1 + r2 + r4 / V0
% dcN2dv = rn2 = r2 + r4 / V0
% dcNO2dv = rno2 = r3 / V0

%creating a vector of volume values from 0 to total volume of the PFR
v_span = linspace(0,V,100);

%calculating the concentrations
[v,c] = ode45(@(v,c)func(v,c,V0,k1,k2,k3,k4),v_span,[cNH30 cO20 0 0 0 0]);

%assigning the correspondent concentrations
cNH3 = c(:,1);
cO2 = c(:,2);
cNO = c(:,3);
cH2O = c(:,4);
cN2 = c(:,5);
cNO2 = c(:,6);

%calculating the molar flows
FNH3 = cNH3 .* V0;
FO2 = cO2 .* V0;
FNO = cNO .* V0;
FH2O = cH2O .* V0;
FN2 = cN2 .* V0;
FNO2 = cNO2 .* V0;

rk = rank(StoichiMat); % rank = 3

StoichiMat2 = StoichiMat;
StoichiMat2(4,:) = []; % removing reaction 4

rk2 = rank(StoichiMat2); % rank is still 3

%calculating the new concentrations
[v,c] = ode45(@(v,c)func2(v,c,V0,k1,k2,k3),v_span,[cNH30 cO20 0 0 0 0]);

%assigning the correspondent concentrations
cNH3_2 = c(:,1);
cO2_2 = c(:,2);
cNO_2 = c(:,3);
cH2O_2 = c(:,4);
cN2_2 = c(:,5);
cNO2_2 = c(:,6);

%calculating the molar flows
FNH3 = cNH3_2 .* V0;
FO2 = cO2_2 .* V0;
FNO = cNO_2 .* V0;
FH2O = cH2O_2 .* V0;
FN2 = cN2_2 .* V0;
FNO2 = cNO2_2 .* V0;

hold on
plot(v_span,cNH3,"Color","#0072BD");
plot(v_span,cO2,"Color","#D95319");
plot(v_span,cNO,"Color","#EDB120");
plot(v_span,cH2O,"Color","#7E2F8E");
plot(v_span,cN2,"Color","#77AC30");
plot(v_span,cNO2,"Color","#f74abb");
plot(v_span,cNH3_2,"Color", "#0072BD", "LineStyle", "--");
plot(v_span,cO2_2,"Color","#D95319", "LineStyle", "--");
plot(v_span,cNO_2,"Color","#EDB120", "LineStyle", "--");
plot(v_span,cH2O_2,"Color","#7E2F8E", "LineStyle", "--");
plot(v_span,cN2_2,"Color","#77AC30", "LineStyle", "--");
plot(v_span,cNO2_2,"Color","#f74abb", "LineStyle", "--");
xlabel("Volume of the reactor (L)");
ylabel("Concentration of component (mol/L)");
legend("NH3","O2","NO","H20","N2","NO2","NH3_2","O2_2","NO_2","H20_2","N2_2","NO2_2");
title("Concentration along the volume of a PFR");
hold off

%with the above plot it can be notted that by removing the dependent
%reaction the new values obtained for concentrations of each component are
%quite diferent, the is a larger concentration of ammonia on the exit of
%the PFR by removing the 4th equation from the simulation for example.

function f = func(v,c,V0,k1,k2,k3,k4)
dcnh3dv = (-k1*c(1)*(c(2)^2) - k2*c(1)*c(2) - k4*c(3)*(c(1)^(2/3)))/V0;
dco2dv = (-k1*c(1)*(c(2)^2) - k2*c(1)*c(2) - k3*(c(3)^2)*c(2))/V0;
dcnodv = (k1*c(1)*(c(2)^2) - k3*(c(3)^2)*c(2) - k4*c(3)*(c(1)^(2/3)))/V0;
dch2odv = (k1*c(1)*(c(2)^2) + k2*c(1)*c(2) + k4*c(3)*(c(1)^(2/3)))/V0;
dcn2dv = (k2*c(1)*c(2) + k4*c(3)*(c(1)^(2/3)))/V0;
dcno2dv = (k3*(c(3)^2)*c(2))/V0;
f = [dcnh3dv;dco2dv;dcnodv;dch2odv;dcn2dv;dcno2dv];
end


function u = func2(v,c,V0,k1,k2,k3)
dcnh3dv = (-k1*c(1)*(c(2)^2) - k2*c(1)*c(2))/V0;
dco2dv = (-k1*c(1)*(c(2)^2) - k2*c(1)*c(2) - k3*(c(3)^2)*c(2))/V0;
dcnodv = (k1*c(1)*(c(2)^2) - k3*(c(3)^2)*c(2))/V0;
dch2odv = (k1*c(1)*(c(2)^2) + k2*c(1)*c(2))/V0;
dcn2dv = (k2*c(1)*c(2))/V0;
dcno2dv = (k3*(c(3)^2)*c(2))/V0;
u = [dcnh3dv;dco2dv;dcnodv;dch2odv;dcn2dv;dcno2dv];
end
