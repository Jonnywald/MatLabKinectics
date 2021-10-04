%Problem set 5
%problem 1

%RAW 6.2

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

Ca0 = 2.0; % mol/L
Cb0 = 2.0; % mol/L
Cc0 = 0.0; % mol/L

km = 0.01725; % L/mol.min
dHr = -10.0; % kcal/mol A @ 27C
Cpa = 20.0; % cal/mol.K
Cpb = 20.0; % cal/mol.K
Cpc = 40.0; % cal/mol.K

Tm = 300.0; %K
EaR = 2660.0; %K

VR = 1200; % L

Tr = 27.0; %C

Q = -41.7 * 1000; %cal/min

% Enrg balance
% 0 = Q - (dHr*dCp(T-Tr))*z Batch
% 0 = Q - dHr*dCp*T*z + dHr*dCp*Tr*z
% dHr*dCp*T*z = Q + dHr*dCp*Tr*z
% T = Q/dHr*dCp*z + Tr

%material Balance
% dCadt = - kCaCb
% dCadt = - kCa²

%Unit conversions
Tr = Tr + 273.15; %converting C to K
dHr = dHr * 1000; %converting kcal/MolA to cal/MolA

dCp = -Cpa -Cpb + Cpc;

t_span = linspace(0,20,100);
[t x] = ode45(@(t,x)fun(t,x,EaR,Tm,Q,dHr,dCp,Tr,km),t_span,[Ca0]);

Ca = x(:,1);

plot(t,Ca);

function f=fun(t,x,EaR,Tm,Q,dHr,dCp,Tr,km)
    
    Ca = x(1);
    z = 2 - Ca;
    T = Q/(dHr*dCp*z) + Tr;
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = - k*Ca^2;
    
    dCadt = ra;
    
    f = [dCadt];

end



