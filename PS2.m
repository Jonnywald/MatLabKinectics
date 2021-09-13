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

Ka1T0 = exp((-dGfMCP)/(R*T0)); %calc of the Ka1 at 600K
Ka2T0 = exp((-dGf2_MP)/(R*T0));%calc of the Ka2 at 600K

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
% 130.1454 K

%composition:
zeta = Solution(1);
yMCP = (zeta)/(1+zeta);     % yMCP = 0.324621838106545
y2_MP = (zeta)/(1+zeta);    % y2_MP = 0.324621838106545
yH2 = (zeta)/(1+zeta);      % yH2 = 0.324621838106545
yHex = (1-2*zeta)/(1+zeta); % yHex = 0.026134485680364

%checking results
yiT = yMCP + y2_MP + yH2 + yHex; %yT = 1


% Problem 2)
% RAW 3.7

% Problem 3)
% RAW 4.2

% Problem 4)



function u = eqn(z)
    u(1) = ((z(1).^2)/(1-z(1)-2*(z(1).^2)))-exp(((31.69)/(0.082057366080960))*((1/z(2))-(1/600))+log(0.394778242347233));
    u(2) = (z(1)/(1-2*z(1)))-exp(((46.10)/(0.082057366080960))*((1/z(2))-(1/600))+log(0.422746304052666));
end
