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


%the clausius-clapeyron equation
%               dHvap
% Ln Pc = c - --------
%                RT

c = 7.53; % clausius-clapeyron constant
T0 = 298; % inital temperature of 298 K
dH0 = -10; % kcal/mol
dHvap = 5; % kcal/mol
K = 8; % K of the reaction
R = 1.98720425864083; % Gas constant kcal * K-1 * mol-1
P = 1.0; %atm

%a) over what temperature range does the reactor contain a liquid phase?
lnPc = c - (dHvap)/(R*T0);
Pc = exp(lnPc);

% Problem 3)
% RAW 4.2

% A -> 2B

% Species  | Feed | Exit (Ni)      |           yi          |     
%----------------------------------------------------------|
% [A]      |  1   |   1 - zeta     | (1 - zeta)/(1 + zeta) |
% [B]      |  0   |    2*zeta      |  (2*zeta)/(1 + zeta)  |
%----------|------|----------------|-----------------------|
% Total    |  1   |  1 + zeta(NiT) |           1           |


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

function u = eqn(z)
    u(1) = ((z(1).^2)/(1-z(1)-2*(z(1).^2)))-exp(((-14.410000000000000)/(0.082057366080960))*((1/z(2))-(1/600))+log(0.945868124008406));
    u(2) = (z(1)/(1-2*z(1)))-exp(((1.530000000000001)/(0.082057366080960))*((1/z(2))-(1/600))+log(1.012878144368649));
end
