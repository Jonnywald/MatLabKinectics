%midterm exam 2
%guilherme Bertola

%A -> B + 0.5 C

% B <-> 2D

%experimental data: time (min) per pressure (torr)
data = [0 268.7; 20 293.0; 40 302.2; 60 311.0; 80 318.9;...
    100 325.9; 120 332.3; 140 338.8; 160 344.4];

%separating the values
t = data(:,1);
P = data(:,2);

%equilibrium constant
Kp = 97.5; %torr

para = 0.0058; %guess
k = nlinfit(t,P,@func,para);
k = real(k); %8.572e-04 min-1
newP = func(k,t);


hold on;
plot(t,P,"k*");
plot(t,newP,"k-");
xlabel("time (min)");
ylabel("pressure (torr)");
title("Pressure v Time");
legend("Experimental","Predicted");
hold off;

function p=func(k,t)
     %initial pressure of N2O5
     Pa0 = 268.7; %torr
     %equilibrium constant
     Kp = 97.5; %torr
     Pa = Pa0 * exp(-k*t);
     Pb = Pa0 - Pa;
     Pc = (Pa0 - Pa)/2;
     Pd = sqrt(Pb*Kp);

     p = Pa + Pb + Pc + Pd;
end