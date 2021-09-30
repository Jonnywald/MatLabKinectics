% Problem set 2
% Guilherme Bertola


% Problem 1)
% RAW 3.4

% Hex <> MCP + H2  (1)
% Hex <> 2_MP    (2)

%where: Hex = n-hexane
%       MCP = Methylcyclopentane
%       2_MP = 2-methylpentane

% Species
% [Hex]
% [MCP]
% [2_MP]
% [H2]

%Heats of reaction (Kcal/mol) at 600K
dHfHex = -46.10;
dHfMCP = -31.69;
dHf2_MP = -47.63;
dHfH2 = 0;

%Gibbs free energy (Kcal/mol) at 600K
dGfHex = 43.02;
dGfMCP = 45.76;
dGf2_MP = 42.39;
dGfH2 = 0;

% Stoichiometric Matrix
StoicMat = [-1 1 0 1;
            -1 0 1 0];
        
T0 = 600; % temperature of the standart state of the System is 600 K
Pt = 1.0; % total pressure of the system in atm
V = 100; % Volume of 100 L
R = 0.082057366080960; % Gas constant in L*atm*K-1*mol-1

% find T where NiMCP = Ni2_MP
% T where zeta1 = zeta2

% Species  | Feed | Exit (Ni)       |           yi          |     
%-----------------------------------------------------------|
% [Hex]    |  1   | 1 -zeta1 -zeta2 | (1 -zeta1 -zeta2)/NiT |
% [MCP]    |  0   |      +zeta1     |   +zeta1/NiT          |
% [2_MP]   |  0   |      +zeta2     |   +zeta2/NiT          |
% [H2]     |  0   |      +zeta1     |   +zeta1/NiT          |
%----------|------|-----------------|-----------------------|
% Total    |  1   | 1 + zeta1 (NiT) |           1           |

dG1 = dGfMCP + dGfH2 - dGfHex;
dG2 = dGf2_MP - dGfHex;

dH1 = dHfMCP + dHfH2 - dHfHex;
dH2 = dHf2_MP - dHfHex;

Ka1T0 = exp((-dG1)/(R*T0)); %calc of the Ka1 at 600K
Ka2T0 = exp((-dG2)/(R*T0));%calc of the Ka2 at 600K

% Ka1 = aMCP * aH2
%       ----------
%           aHex

% Ka2 =    a2_MP
%       ----------
%          aHex

% Ka1 = yMCP * yH2  
%       ----------  *  Pt/1atm
%           yHex

% Ka2 =    y2_MP
%       ---------- *  Pt/1atm
%          yHex

% Ka1 =           zeta1 * zeta1  
%       -------------------------------  *  Pt/1atm
%       (1 - zeta1 -zeta2) * (1+ zeta1)

% Ka2 =        zeta2
%       ------------------ *  Pt/1atm
%       (1 - zeta1 -zeta2)

%If zeta1 = zeta2

% Ka1 =               (zeta1)^2  
%       ----------------------------------  *  Pt/1atm
%           1 + zeta1 - 2zeta1 -(2zeta1)^2

% Ka2 =        zeta1
%       ------------------ *  1atm/1atm
%          (1 - 2zeta1)

% using van hoff

%    0(1) = ((z.^2)/(1-z-2*(z.^2)))-exp(((-dHfMCP)/(R))*((1/T)-(1/T0))+log(Ka1T0));
%    0(2) = (z/(1-2*z))-exp(((-dHf2_MP)/(R))*((1/T)-(1/T0))+log(Ka2T0));

x0 = [0.4,100]; %initial guess for the fsolve
Solution = fsolve(@eqn,x0); %solution for the system

% the temperature where equal amounts of MCP and 2_MP are formed is
% 120.710626318075 K

%composition:
zeta = Solution(1);
yMCP = (zeta)/(1+zeta);     % yMCP = 0.258219006548342
y2_MP = (zeta)/(1+zeta);    % y2_MP = 0.258219006548342
yH2 = (zeta)/(1+zeta);      % yH2 = 0.258219006548342
yHex = (1-2*zeta)/(1+zeta); % yHex = 0.225342980354973

%checking results
yiT = yMCP + y2_MP + yH2 + yHex; %yT = 1


% Problem 2)
% RAW 3.7

% gas phase reaction
% A + B <> C

% Species  | Feed | Exit (Ni)      |           yi          |     
%----------------------------------------------------------|
% [A]      |  1   |   1 - zeta     | (1 - zeta)/(2 - zeta) |
% [B]      |  1   |   1 - zeta     | (1 - zeta)/(2 - zeta) |
% [C]      |  0   |     zeta       |   (zeta)/(2 - zeta)   |
%----------|------|----------------|-----------------------|
% Total    |  2   | 2 - zeta (NiT) |           1           |

%         yC
% Ka = --------- * Pt/1atm
%       yA * yB

%                zeta
% Ka = ------------------------ * (2 - zeta) * Pt/1atm
%       (1 - zeta) * (1 - zeta)

%the clausius-clapeyron equation
%               dHvap
% Ln Pc = c - --------
%                RT

c = 7.53; % clausius-clapeyron constant
T0 = 298; % inital temperature of 298 K
dH0 = -10; % kcal/mol
dHvap = 5; % kcal/mol
K = 8; % K of the reaction
R = 1.98720425864083 / 1000; % Gas constant kcal * K-1 * mol-1
Pt = 1.0; %atm
V = 100; %L
%a) over what temperature range does the reactor contain a liquid phase?
lnPc = c - (dHvap)/(R*T0);
Pc = exp(lnPc);

% Vapor pressure of C is 0.4012 atm

% T = Pc*V/N*R
% T = (Pc*V)/((2-z)*R)
% 0 = (Pc*V)/((2-z)*R) - T
% van hoff

% Ln (Ka/8) = (-dH0/R)*(1/T - 1/T0)
% Ln Ka - Ln 8 = (-dH0/R)*(1/(T)-1/T0)
% Ka = exp((-dH0/R)*(1/T-1/T0) + ln(8))

% with the activities:
% Ka = (z/((1-z)*(1-z)))*(2-z)

% 0 =  exp((-dH0/R)*(1/(T)-1/T0) + log(8)) - (z/((1-z)*(1-z)))*(2-z)

x0 = [0.6,500];
val = fsolve(@(z)prob2(z,dH0,R,T0,Pc,V),x0);
z = val;
%calculatin the Temp where Pt = Pc
T = (Pc*V)/((2-z(1))*R);

% T = 2.026687810582085e+03 K

%b)
dH0 = 10; %Kcal/mol
x0 = [0.6,500];
val = fsolve(@(z)prob2(z,dH0,R,T0,Pc,V),x0);
z = val;
T = (Pc*V)/((2-z(1))*R);

% T = 1.604381156135477e+02 K

% Problem 3)
% RAW 4.2

% A -> 2B

% Species  | Feed | Exit (Ni)      |           yi          |     
%----------------------------------------------------------|
% [A]      |  1   |   1 - zeta     | (1 - zeta)/(1 + zeta) |
% [B]      |  0   |    2*zeta      |  (2*zeta)/(1 + zeta)  |
%----------|------|----------------|-----------------------|
% Total    |  1   |  1 + zeta(NiT) |           1           |
 
k = 0.35; % min-1
Pt = 1.0; % atm
ca0 = 1;
t = 5; %min

%Lncaca0 = -kt
%ca = exp(-kt + ln(ca0))

ca = exp(-k*t + log(ca0));

%ca = 0.173773943450445

% Problem 4)

% C(s) + 2H2O(g) <> CO2(g) + 2H2(g) (1)

% C(s) + H2O(g)  <> CO(g) + H2(g)   (2)

% C(s) + 2H2(g)  <> CH4             (3)

% Species  | Feed | Exit (Ni)       |           yi          |     
%-----------------------------------------------------------|
% [H2O]    |  1   |  1 - 2*z1 - z2  |  (1 - 2*z1 - z2)/NiT  |
% [CO2]    |  0   |      +z1        |      +z1/NiT          |
% [CO]     |  0   |      +z2        |      +z2/NiT          |
% [H2]     |  0   | 2*z1 + z2 -2*z3 |(2*z1 + z2 -2*z3)/NiT  |
% [CH4]    |  0   |      +z3        |      +z3/NiT          |
%----------|------|-----------------|-----------------------|
% Total    |  1   | 1+z1+z2+z3 (NiT)|           1           |

%       yCO2 * yH2^2    1bar
% Ka1 = ------------- * ----
%           yH2O^2      1bar

%         yCO * yH2     1bar
% Ka2 = ------------- * ----
%           yH2O        1bar

%         yCH4    1bar
% Ka3 = ------- * ----
%        yH2^2    1bar

%       z1 * (2*z1 + z2 -2*z3)^2          1
% Ka1 = ------------------------- * -------------
%           (1 - 2*z1 - z2)^2       (1+z1+z2+z3)

%         z2 * (2*z1 + z2 -2*z3)           1
% Ka2 = -------------------------- * -------------
%           (1 - 2*z1 - z2)          (1+z1+z2+z3)

%                z3    
% Ka3 = --------------------- * (1+z1+z2+z3)
%        (2*z1 + z2 -2*z3)^2

T0 = 298.15; %standart state temp in K
R = 8.31446261815324 / 1000; % Gas constant KJ . K-1 . mol-1

% heat of formation of each component at standart state:
dH_H2O = -241.8; %KJ/mol
dH_CO2 = -393.5; %KJ/mol
dH_CO = -110.5; %KJ/mol
dH_H2 = 0; %KJ/mol
dH_CH4 = -74.5; %KJ/mol

% gibbs free energy of each component at standart state:
dGH2O = -228.6; %KJ/mol
dGCO2 = -394.4; %KJ/mol
dGCO = -137.2; %KJ/mol
dGH2 = 0; %KJ/mol
dGCH4 = -50.5; %KJ/mol

%heat capacities of the components:
aH2O = 32.218;
aCO2 = 22.243;
aCO = 28.142;
aH2 = 29.088;
aCH4 = 19.875;
bH2O = 0.192 / 100;
bCO2 = 5.977 / 100;
bCO = 0.167 / 100;
bH2 = -0.192 / 100;
bCH4 = 5.021 / 100;
cH2O = 1.055 / 100000;
cCO2 = -3.499 / 100000;
cCO = 0.537 / 100000;
cH2 = 0.4 / 100000;
cCH4 = 1.268 / 100000;
dH2O = -3.593 / 1000000000;
dCO2 = 7.464 / 1000000000;
dCO = -2.221 / 1000000000;
dH2 = -0.870 / 1000000000;
dCH4 = -11.004 / 1000000000;

% calc of the delta "a" of each reaction
da1 = -2*aH2O + aCO2 + 2*aH2;
da2 = -aH2O + aCO + aH2;
da3 = -2*aH2 + aCH4;

% calc of the delta "b" of each reaction
db1 = -2*bH2O + bCO2 + 2*bH2;
db2 = -bH2O + bCO + bH2;
db3 = -2*bH2 + bCH4;

% calc of the delta "c" of each reaction
dc1 = -2*cH2O + cCO2 + 2*cH2;
dc2 = -cH2O + cCO + cH2;
dc3 = -2*cH2 + cCH4;

% calc of the delta "d" of each reaction
dd1 = -2*dH2O + dCO2 + 2*dH2;
dd2 = -dH2O + dCO + dH2;
dd3 = -2*dH2 + dCH4;

% calc of the delta H of the standart state for each reaction
dHt0_1 = -2*dH_H2O + dH_CO2 + 2*dH_H2;
dHt0_2 = -dH_H2O + dH_CO + dH_H2;
dHt0_3 = -2*dH_H2 + dH_CH4;

% calc of the delta G of the standart state for each reaction
dGt0_1 = -2*dGH2O + dGCO2 + 2*dGH2;
dGt0_2 = -dGH2O + dGCO + dGH2;
dGt0_3 = -2*dGH2 + dGCH4;

%calc of Kt0 of each reaction:
Kt0_1 = exp((-dGt0_1)/(R*T0));
Kt0_2 = exp((-dGt0_2)/(R*T0));
Kt0_3 = exp((-dGt0_3)/(R*T0));

%settng up the temperature range of 600 K to 1600K
t_range = linspace(600,1600,1000); 

%initial Guess for fsolve
x0 = [0.05,0.05,0.05];
i = 1;

%Preallocation of matrix to save time
molH2O = zeros(1,1000);
molCO2 = zeros(1,1000);
molCO = zeros(1,1000);
molH2 = zeros(1,1000);
molCH4 = zeros(1,1000);

for t = t_range
    
    %item 1 of van hoff for the reaction 1
    item1_1 = (da1/R)*log(t/T0);
    %item 2 of van hoff for the reaction 1
    item2_1 = (db1/(2*R))*(t-T0);
    %item 3 of van hoff for the reaction 1
    item3_1 = (dc1/(6*R))*((t.^2)-(T0.^2));
    %item 4 of van hoff for the reaction 1
    item4_1 = (dd1/(12*R))*((t.^3)-(T0.^3));
    %item 5 of van hoff for the reaction 1
    item5_1 = (1/R)*(-dHt0_1+da1*T0+(db1/2)*(T0.^2)+(dc1/3)*(T0.^3)+(dd1/4)*(T0.^4))*((1/t)-(1/T0));
    
    %item 1 of van hoff for the reaction 2
    item1_2 = (da2/R)*log(t/T0);
    %item 2 of van hoff for the reaction 2
    item2_2 = (db2/(2*R))*(t-T0);
    %item 3 of van hoff for the reaction 2
    item3_2 = (dc2/(6*R))*((t.^2)-(T0.^2));
    %item 4 of van hoff for the reaction 2
    item4_2 = (dd2/(12*R))*((t.^3)-(T0.^3));
    %item 5 of van hoff for the reaction 2
    item5_2 = (1/R)*(-dHt0_2+da2*T0+(db2/2)*(T0.^2)+(dc2/3)*(T0.^3)+(dd2/4)*(T0.^4))*((1/t)-(1/T0));
    
    %item 1 of van hoff for the reaction 3
    item1_3 = (da3/R)*log(t/T0);
    %item 2 of van hoff for the reaction 3
    item2_3 = (db3/(2*R))*(t-T0);
    %item 3 of van hoff for the reaction 3
    item3_3 = (dc3/(6*R))*((t.^2)-(T0.^2));
    %item 4 of van hoff for the reaction 3
    item4_3 = (dd3/(12*R))*((t.^3)-(T0.^3));
    %item 5 of van hoff for the reaction 3
    item5_3 = (1/R)*(-dHt0_3+da3*T0+(db3/2)*(T0.^2)+(dc3/3)*(T0.^3)+(dd3/4)*(T0.^4))*((1/t)-(1/T0));
    
    %calc of the van hoff eqn for each reaction
    lnktkt0_1 = item1_1 + item2_1 + item3_1 + item4_1 + item5_1;
    lnktkt0_2 = item1_2 + item2_2 + item3_2 + item4_2 + item5_2;
    lnktkt0_3 = item1_3 + item2_3 + item3_3 + item4_3 + item5_3;
    
    %using the aproximation since there may be a bug in the code
    lnktkt0_1 = (-dHt0_1/R)*(1/t-1/T0);
    lnktkt0_2 = (-dHt0_2/R)*(1/t-1/T0);
    lnktkt0_3 = (-dHt0_3/R)*(1/t-1/T0);
    
    % calculating the K for each reaction
    Kt1 = exp(lnktkt0_1+log(Kt0_1));
    Kt2 = exp(lnktkt0_2+log(Kt0_2));
    Kt3 = exp(lnktkt0_3+log(Kt0_3));
       
    %solving the system
    values = fsolve(@(z)prob4(z,Kt1,Kt2,Kt3),x0);
    
    %assinging the extent of reactions
    z1 = values(1);
    z2 = values(2);
    z3 = values(3);
    
    %calculating the composition:
    NiT = 1 + z1 + z2 + z3;
    yH2O = (1 - 2*z1 - z2)/NiT;
    yCO2 = z1/NiT;
    yCO = z2/NiT;
    yH2 = (2*z1 + z2 -2*z3)/NiT;
    yCH4 = z3/NiT;
    
    yT = yH2O + yCO2 + yCO + yH2 + yCH4;
    
    %creating a vector with the molar fraction of each component:
    molH2O(i) = yH2O;
    molCO2(i) = yCO2;
    molCO(i) = yCO;
    molH2(i) = yH2;
    molCH4(i) = yCH4;
    
    i = i + 1;
end

%plot for H2O
plot(t_range,molH2O,"k");
xlabel("Temperature (K)");
ylabel("Molar composition");
title("yH2O over temperature");

%plot for CO2
plot(t_range,molCO2,"k");
xlabel("Temperature (K)");
ylabel("Molar composition");
title("yCO2 over temperature");

%plot for CO
plot(t_range,molCO,"k");
xlabel("Temperature (K)");
ylabel("Molar composition");
title("yCO over temperature");

%plot for H2
plot(t_range,molH2,"k");
xlabel("Temperature (K)");
ylabel("Molar composition");
title("yH2 over temperature");

%plot for CH4
plot(t_range,molCH4,"k");
xlabel("Temperature (K)");
ylabel("Molar composition");
title("yCH4 over temperature");


function F = prob4(z,Kt1,Kt2,Kt3)
    F(1) = ((z(1)*((2*z(1)+z(2)-2*z(3)).^2))/((1-2*z(1)-z(2)).^2))*(1/(1+z(1)+z(2)+z(3))) - Kt1;
    F(2) = ((z(2)*(2*z(1)+z(2)-2*z(3)))/(1-2*z(1)-z(2)))*(1/(1+z(1)+z(2)+z(3)))- Kt2;
    F(3) = (z(3)/((2*z(1)+z(2)-2*z(3)).^2))*(1+z(1)+z(2)+z(3))- Kt3;
end
function w = prob2(z,dH0,R,T0,Pc,V)
    w(1) = exp((-dH0/R)*(1/(z(2))-1/T0) + log(8)) - (z(1)/((1-z(1))*(1-z(1))))*(2-z(1));
    w(2) = (Pc*V)/((2-z(1))*R) - z(2);
end
function u = eqn(z)
    u(1) = ((z(1).^2)/(1-z(1)-2*(z(1).^2)))-exp(((-14.410000000000000)/(0.082057366080960))*((1/z(2))-(1/600))+log(0.945868124008406));
    u(2) = (z(1)/(1-2*z(1)))-exp(((1.530000000000001)/(0.082057366080960))*((1/z(2))-(1/600))+log(1.012878144368649));
end
