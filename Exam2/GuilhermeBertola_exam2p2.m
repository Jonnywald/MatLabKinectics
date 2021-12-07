%midterm exam 2
%guilherme Bertola

%A -> B + 0.5 C

% B <-> 2D

%experimental data: time (min) per pressure (torr)
data = [0 268.7; 20 293.0; 40 302.2; 60 311.0; 80 318.9;...
    100 325.9; 120 332.3; 140 338.8; 160 344.4];
%equilibrium pressure
Eq_pressure = 473.0;

%separating the values
t = data(:,1);
P = data(:,2);

%temperature
T = 298; %K

%initial pressure of N2O5
Pa0 = 268.7; %torr



%equilibrium constant
K = 97.5; %torr

k = nlinfit(t,P,@func,K);

newP = func(k,t);

hold on;
plot(t,P,"k*");
plot(t,newP,"k-");
xlabel("time (min)");
ylabel("pressure (torr)");
hold off;

%-ra = KPa

% Pa = 
function p=func(k,t)
    p = 1.5*268.7 - 0.5*268.7*exp(-k*t);
end