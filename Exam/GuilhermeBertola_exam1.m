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
%Ft0 = Ft

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
Ua = 16.0; %J/(min*dm3*C) = J/(min*L*C)

%Start of Calculations

%Initial concentrations
Ca0 = Ct0 * ya0;
Cb0 = Ct0 * yb0;
Cc0 = 0.0;
%Initial mol flow rates
Fa0 = Ft0 * ya0;
Fb0 = Ft0 * yb0;

%calc initial volumetric flow
V0 = Ft0/Ct0;

%defining a Volume to be used in the calculations
Vr = 100; %L
%Volume of the reactor (Vr)
v_span = linspace(0,Vr,1000);
%initial conditions of the ODE
y0 = [Ca0 Cb0 Cc0 T0];
%calling ODE45 to solve the exercise
[v,x] = ode45(@(v,x)PFR(v,x,V0,T0,Ua,Ta,dHr1,dHr2,Cpa,Cpb,Cpc),v_span,y0);

%assigning the values
Ca = x(:,1);
Cb = x(:,2);
Cc = x(:,3);
T = x(:,4);

%A)
%plotting the graphs
subplot(2,2,1);
plot(v,Ca);
xlabel("Volume (L)");
ylabel("Concentration of A (mol/L)");
title("Ca x Volume");
subplot(2,2,2);
plot(v,Cb);
xlabel("Volume (L)");
ylabel("Concentration of B (mol/L)");
title("Cb x Volume");
subplot(2,2,3);
plot(v,Cc);
xlabel("Volume (L)");
ylabel("Concentration of C (mol/L)");
title("Cc x Volume");
subplot(2,2,4);
plot(v,T);
xlabel("Volume (L)");
ylabel("Temperature (K)");
title("T x Volume");
sgtitle("Plots for part ""a"" ");

%B)
%find the lowest concentration of species A in the reactor
%assuming our reactor have a volume of Vr = 100L
LowestA = min(Ca);

%C)
%find the Highest concentration of species B in the reactor
%assuming our reactor have a volume of Vr = 100L
HighestB = max(Cb);

%D)
%find the Highest concentration of species B in the reactor
%assuming our reactor have a volume of Vr = 100L
HighestA = max(Ca);

%E)
%repeat for a pure feed of A
Ca0 = 2.0;
Cb0 = 0.0;

%initial conditions of the ODE
y0 = [Ca0 Cb0 Cc0 T0];
%calling ODE45 to solve the exercise
[v,x] = ode45(@(v,x)PFR(v,x,V0,T0,Ua,Ta,dHr1,dHr2,Cpa,Cpb,Cpc),v_span,y0);

%assigning the values
Ca = x(:,1);
Cb = x(:,2);
Cc = x(:,3);
T = x(:,4);

%A - E)
%plotting the graphs
subplot(2,2,1);
plot(v,Ca);
xlabel("Volume (L)");
ylabel("Concentration of A (mol/L)");
title("Ca x Volume");
subplot(2,2,2);
plot(v,Cb);
xlabel("Volume (L)");
ylabel("Concentration of B (mol/L)");
title("Cb x Volume");
subplot(2,2,3);
plot(v,Cc);
xlabel("Volume (L)");
ylabel("Concentration of C (mol/L)");
title("Cc x Volume");
subplot(2,2,4);
plot(v,T);
xlabel("Volume (L)");
ylabel("Temperature (K)");
title("T x Volume");

sgtitle("Plots for part ""e"" ");
%B - E)
%find the lowest concentration of species A in the reactor
%assuming our reactor have a volume of Vr = 100L
LowestA_E = min(Ca);

%C - E)
%find the Highest concentration of species B in the reactor
%assuming our reactor have a volume of Vr = 100L
HighestB_E = max(Cb);

%D - E)
%find the Highest concentration of species B in the reactor
%assuming our reactor have a volume of Vr = 100L
HighestA_E = max(Ca);

%the function for the exercise
function f=PFR(v,x,V0,T0,Ua,Ta,dHr1,dHr2,Cpa,Cpb,Cpc)
    %setting up the variables
    Ca = x(1);
    Cb = x(2);
    Cc = x(3);
    T = x(4);
    %volumetric flow
    V = V0*(T/T0);
    %reactions K
    kc = 10.0*exp(4.8*(430/T-1.5));
    k1 = 0.5*exp(2*(1-320/T));
    k2 = k1/kc;
    k3 = 0.005*exp(4.6*(1-460/T));
    %heat exchange
    %since dT in K = dT in C should be fine to use the values in K even if
    %Ua has a unit of C
    Q = Ua*(T-Ta);
    %reaction rates
    r1 = k1*Ca - k2*Cb;
    r2 = k3*Ca;
    %components rates
    ra = -r1-r2;
    rb = r1;
    rc = r2;
    % the mol balance and energy balance (ODES)
    dCadv = ra/V;
    dCbdv = rb/V;
    dCcdv = rc/V;
    dTdv = ((-r1)*(-dHr1)+(-r2)*(-dHr2)-Q)/((Ca*Cpa+Cb*Cpb+Cc*Cpc)*V);
    
    f = [dCadv;dCbdv;dCcdv;dTdv];
end

