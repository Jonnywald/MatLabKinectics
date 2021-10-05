%Problem set 5
% problem 4

% Handout P 8-9 (a-e)

% 2 VinylAcetylene -> Stylene
% 2A -> B

% Species | MolFeed Fi0| Exit (Fi) |  Exit    |     Ci      |
%---------|------------|-----------|----------|-------------|
% A       |    Fa0     |    Fa     |  Fa0 -2z | (Fa0 -2z)/V |
% B       |    0       |    Fb     |   + z    |    z/V      |
%---------|------------|-----------|----------|-------------|
% Total   |Ff = Fa0+Fb0|    FT     |   Fa0-z  |   Ff/V      |

Ca0 = 1.0; % mol/dm3
Fa0 = 5.0; % mol/s

V = Fa0/Ca0;

% dHr = -231-0.012(T - 298)  %Kj/Mol

Cpa = 0.1222; %Kj/mol*K

%K = 1.48e11 * exp(-19.124/T) %dm3/mol*s

T0 = 675; %K

Ua = 5.0; %Kj/s*dm3/K
R = 8.31446261815324 / 1000;

Ta = 700; %K

%a) determine the Xa in a Vr = 10 dm3 PFR with T0 = 675
Vr = 10; %dm3
v_span = linspace(0,Vr,200);
y0 = [0 T0];
[v x] = ode45(@(v,x)part_a(v,x,Fa0,Cpa,Ta,Ua,Ca0),v_span,y0);
Xa = x(:,1);
temp = x(:,2);


subplot(1,2,1);
plot(v,Xa);
xlabel("Volume (dm3)");
ylabel("Conversion Xa");
title("Conversion x Volume");
subplot(1,2,2);
plot(v,temp);
xlabel("Volume (dm3)");
ylabel("temperature");
title("temperature x Volume");

function f=part_a(v,x,Fa0,Cpa,Ta,Ua,Ca0)
    Xa = x(1);
    T = x(2);
    Ca = Ca0  - Xa*Ca0;
    dHr = -231-0.012*(T - 298);
    k = 1.48e11 * exp(-19.124/T);
    ra = k*Ca^2;
    dXadv = ra/Fa0;
    dTdv = (Ua*(Ta-T)-ra*dHr)/(Fa0*Cpa);
    f = [dXadv;dTdv];
end
