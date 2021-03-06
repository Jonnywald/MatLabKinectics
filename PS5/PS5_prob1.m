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

t_span = linspace(0,15,100);
[t x] = ode45(@(t,x)fun(t,x,EaR,Tm,Q,dHr,km,VR,Cpa),t_span,[Ca0 Tr]);

Ca = x(:,1);
Xa = (Ca0 - Ca)./Ca0;
T = x(:,2);

subplot(1,3,1);
plot(t,Ca);
xlabel("time (min)");
ylabel("Concetration Ca (mol/L)");
title("Ca x Time");
subplot(1,3,2);
plot(t,Xa);
xlabel("time (min)");
ylabel("Conversion Xa");
title("Xa x Time");
subplot(1,3,3);
plot(t,T);
xlabel("time (min)");
ylabel("Temperature (K)");
title("Temp x Time");

%a)
tmax = max(T); %2765.18 K
timeTmax = t(find(T==tmax,1)); %6.21 min

%b)
i=1;
for temp=T'
   if temp<=Tr & i~=1
       indexTimeT_return = i;
       break;
   end
   i=i+1;
end

timeT_return = t(indexTimeT_return); %13.18 min

%c)
i=1;
for time=Xa'
   if time>=0.95
       iXa = i;
       break;
   end
   i=i+1;
end
time95 = t(iXa); %3.63 min

function f=fun(t,x,EaR,Tm,Q,dHr,km,Vr,Cpa)
    
    Ca = x(1);
%     z = 2 - Ca;
%     T = Q/(dHr*dCp*z) + Tr;
    T = x(2);
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = - k*Ca^2;
    
    dCadt = ra;
    dTdt = (-ra*Vr*-dHr+Q)/(Ca*Vr*Cpa);
    f = [dCadt;dTdt];

end



