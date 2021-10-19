%Guilherme Bertola
%Take home portion of midterm exam 1

% Reaction

% A <-> B
% A  -> C

% Species | Fi0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       |   1  |    Fa     | 1-z1-z2  |
% B       |   1  |    Fb     |   1+z1   |
% C       |   0  |    Fc     |   z2     |
%---------|------|-----------|----------|
% Total   |   2  |    Ft     |    2     |

%equimolar feed of A and B
ya0 = 0.5;
yb0 = 0.5;

%Total feed rate
Ft0 = 2.0; % mol/min

%useful information
% all Heat capacities are bearly the same at 100J/mol.K
Cpa = 100.0; %J/Mol/K
Cpb = 100.0; %J/Mol/K
Cpc = 100.0; %J/Mol/K

%concentration of total feed
Ct0 = 2.0; %mol/L

% delta H of rxn
dHr1 = -1800; %J/molA
dHr2 = -1100; %J/molA

%K constants
%k1 = 0.5*exp(2*(1-320/T)); % min-1
%k2 = k1/kc
%k3 = 0.005*exp(4.6*(1-460/T)); % min-1
%kc = 10.0*exp(4.8*(430/T-1.5));

%Initial temperature
T0 = 330.0; %K

%Heat exchange data
Ta = 500.0; %K
Ua = 16.0; %J/(min*dm2*C)

%Start of Calculations

%Initial concentrations
Ca0 = Ct0 * ya0;
Cb0 = Ct0 * yb0;

%Initial mol flow rates
Fa0 = Ft0 * ya0;
Fb0 = Ft0 * yb0;

%calc initial volumetric flow
V0 = Ft0/Ct0;



%the function for the exercise

function f=PFR()
    %setting up the variables
    Ca = x(1);
    Cb = x(2);
    Cc = x(3);
    T = x(4);
    
    %reactions K
    k1 = 0.5*exp(2*(1-320/T));
    k2 = k1/kc;
    k3 = 0.005*exp(4.6*(1-460/T));
    kc = 10.0*exp(4.8*(430/T-1.5));
    %heat exchange
    Q = (2/R)*Ua*(T-Ta);
    %reaction rates
    r1 = k1*Ca;
    r2 = k2*Cb;
    r3 = k3*Ca;
    %components rates
    ra = -r1+r2-r3;
    rb = r1-r2;
    rc = r3;
    % the mol balance and energy balance (ODES)
    dCadv = ra/V;
    dCbdv = rb/V;
    dCcdV = rc/V;
    dTdv = 0; %% ra(dhr)-Q / ficpi
    
    f = [dCadv;dCbdv;dCcdv;dTdv];
end

