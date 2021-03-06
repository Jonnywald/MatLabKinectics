%Problem Set 8
%Guilherme Del Rio

%Problem 2
%Integral Method of Rate Data Analysis in a Batch Reactor

%raw data time (s) vs Pressure (mmHg)
dataset = [0 420; 57 584; 85 662; 114 743; 145 815; 182 891; 219 954;...
    261 1013; 299 1054];

time = dataset(:,1);
P = dataset(:,2);
P0 = P(1);

y = log(((3.*P0)-P)./(2.*P0));

figure;
plot(time,y,"k*");
xlabel("time (s)");
ylabel("Ln(3P0-P/2P0)");
title("Integral method of rate data analysis");

%linear regression
k = 0.004808; %s^-1

%non-linear regression
x = time;
y = P;
val0 = [k];
val = nlinfit(x,y,@func,val0);
newy = func(val,x);

figure;
hold on;
plot(x,y,"k*");
plot(x,newy,"k-");
xlabel("Time (s)");
ylabel("Pressure (mmHg)");
title("Non-linear regression");
hold off;

k = val; % 0.00446 s^-1

function p=func(val,t)
    p = 3*420 - 2*420*exp(-val*t);
end