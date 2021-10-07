%PROBLEM SET 6
%PROBLEM 1
%RAW 6.1

% A + B -> C

% Species | Ci0  | Exit (Fi) |  Exit    |     
%---------|------|-----------|----------|
% A       |   2  |    Ca     |  2 - z   |
% B       |   2  |    Cb     |  2 - z   |
% C       |   0  |    Cc     |  0 + z   |
%---------|------|-----------|----------|
% Total   |   4  |    Ct     |  4 - z   |

% r = kCaCb
% Ca = 2 - z
% z = 2 - Ca

% r = kCa^2

Ca0 = 2.0; % mol/L
Cb0 = 2.0; % mol/L
Cc0 = 0.0; % mol/L

Cpa = 20.0; % cal/mol.K
Cpb = 20.0; % cal/mol.K
Cpc = 40.0; % cal/mol.K

dHr = -10.0; % kcal/mol A @ 27C
VR = 1200; % L
Tr = 27.0; %C
Q = -41.7 * 1000; %cal/min
km = 0.01725; %L/mol/min
Tm = 300; %K
EaR = 2660; %K

%Unit conversions
Tr = Tr + 273.15; %converting C to K
dHr = dHr * 1000; %converting kcal/MolA to cal/MolA

%Material Balance
% dCadt = - kCa²

%Enrg. Balance
% dTdt = (-ra*Vr*-dHr)/(Ca*Vr*Cpa)

y0 = [Ca0, Tr];
t_span = linspace(0,10,500);
[t,x] = ode45(@(t,x)part_a(t,x,km,EaR,Tm,VR,Cpa,dHr),t_span,y0);

Ca = x(:,1);
Temp = x(:,2);
subplot(1,2,1);
plot(t,Ca);
subplot(1,2,2);
plot(t,Temp);
function f=part_a(t,x,km,EaR,Tm,VR,Cpa,dHr)
    Ca = x(1);
    T = x(2);
    
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = -k * Ca^2;
    dCadt = ra;
    dTdt = (-ra*VR*-dHr)/(Ca*VR*Cpa);
    f = [dCadt;dTdt];
end

