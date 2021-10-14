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


%A)
y0 = [Ca0, Tr];
t_span = linspace(0,10,500);
[t,x] = ode45(@(t,x)part_a(t,x,km,EaR,Tm,VR,Cpa,dHr),t_span,y0);

Ca = x(:,1);
Temp = x(:,2);

Temp_rise = Temp(length(Temp))-Temp(1); %3266.5 K

%B)
t_span = linspace(0,1000,1001);
[t,x] = ode45(@(t,x)part_b(t,x,km),t_span,Ca0);
Ca = x(:,1);

Xa = (Ca0 - Ca)./Ca0;

for i=1:1001
    if Xa(i)>=0.95
        Time = t(i); %552 min
        break
    end
end

%time to achive 0.95
% 552 Min

%C)
y0 = [Ca0, Tr];
t_span = linspace(0,10,500);
[t,x] = ode45(@(t,x)part_c(t,x,km,EaR,Tm,VR,Cpa,dHr),t_span,y0);

Ca = x(:,1);
Temp = x(:,2);

Xa = (Ca0 - Ca)./Ca0;
for i=1:500
    if Xa(i)>=0.95
        Time = t(i); %3.6 min
        break
    end
end
plot(Temp,Ca);
xlabel("Temperature (K)");
ylabel("Concentration (mol/L)");
title("Temperature x Concentration of A");

%D)
Ua = 0.01 * 1000; % Cal/(min L K)
Ta = 300; %K
y0 = [Ca0, Tr];
t_span = linspace(0,3,500);
[t,x] = ode45(@(t,x)part_d(t,x,km,EaR,Tm,VR,Cpa,dHr,Ta,Ua),t_span,y0);

Ca = x(:,1);
Temp = x(:,2);

plot(Temp,Ca);
xlabel("Temperature (K)");
ylabel("Concentration (mol/L)");
title("Temperature x Concentration of A");

%E)
Ua = -22.7; % Cal/(min L K)
Ta = 300; %K
y0 = [Ca0, Tr];
t_span = linspace(0,3,500);
[t,x] = ode45(@(t,x)part_e(t,x,km,EaR,Tm,VR,Cpa,dHr,Ta,Ua),t_span,y0);

Ca = x(:,1);
Temp = x(:,2);

tmax = max(Temp); %349.75 K

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

function f=part_b(t,x,k)
    Ca = x(1);
    ra = -k * Ca^2;
    dCadt = ra;

    f = [dCadt];
end

function f=part_c(t,x,km,EaR,Tm,VR,Cpa,dHr)
    Ca = x(1);
    T = x(2);
    
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = -k * Ca^2;
    dCadt = ra;
    dTdt = ((-ra)*VR*(-dHr))/(Ca*VR*Cpa);
    f = [dCadt;dTdt];
end

function f=part_d(t,x,km,EaR,Tm,VR,Cpa,dHr,Ta,Ua)
    Ca = x(1);
    T = x(2);
    
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = -k * Ca^2;
    Q = Ua*(T - Ta)*VR;
    dCadt = ra;
    dTdt = (Q + (-ra*VR*-dHr))/(Ca*VR*Cpa);
    f = [dCadt;dTdt];
end
function f=part_e(t,x,km,EaR,Tm,VR,Cpa,dHr,Ta,Ua)
    Ca = x(1);
    T = x(2);
    
    k = km*exp(-EaR*(1/T - 1/Tm));
    ra = -k * Ca^2;
    Q = Ua*(T - Ta)*VR;
    dCadt = ra;
    dTdt = (Q + (-ra*VR*-dHr))/(Ca*VR*Cpa);
    f = [dCadt;dTdt];
end

