% Problem set 8
% Guilherme Bertola

% PROBLEM 1 - DIFFERENTIAL METHOD OF RATE DATA ANALYSIS
% fit data to a polynomial and diffentiate the polynomial to get rate data
% linearize to extract k and n from the rate expression. repeat using
% nonlinear regression.

% time (min) by Concentration (mol/L)
dataset = [0 0.3335; 2.25 0.2965; 4.50 0.2660; 6.33 0.2450; 8.00 0.2255;...
    10.25 0.2050; 12.00 0.1910; 13.50 0.1794; 15.60 0.1632; 17.85 0.15;...
    19.60 0.1429; 27.00 0.1160; 30.00 0.1053; 38.00 0.0830; 41.00 0.0767;...
    45.00 0.0705; 47.00 0.0678; 57.00 0.0553; 63.00 0.0482];

figure;
plot(dataset(:,1),dataset(:,2),"ko");
xlabel("time (min)");
ylabel("Concentration of Br(mol/L)");
title("Raw data time x Concentration (CBr)");

% polyfit
x1 = dataset(:,1);
m = max(x1);
h = 50;
x2 = linspace(0,m,h);
% 2nd order
% n=2;
% p=polyfit(dataset(:,1),dataset(:,2),n);
% fit1=p(1)*x1.^n + p(2)*x1.^(n-1) + p(n+1);
% fit2=p(1)*x2.^n + p(2)*x2.^(n-1) + p(n+1);
% figure;
% hold on;
% plot(dataset(:,1),dataset(:,2),"ko");
% plot(x1,fit1,"k-");
% xlabel("time (min)");
% ylabel("Concentration of Br(mol/L)");
% title("second order fit");
% hold off;

% 3rd order
% n=3;
% p=polyfit(dataset(:,1),dataset(:,2),n);
% fit1=p(1)*x1.^n + p(2)*x1.^(n-1) + p(3)*x1.^(n-2) + p(n+1);
% fit2=p(1)*x2.^n + p(2)*x2.^(n-1) + p(3)*x2.^(n-2) + p(n+1);
% figure;
% hold on;
% plot(dataset(:,1),dataset(:,2),"ko");
% plot(x1,fit1,"k-");
% xlabel("time (min)");
% ylabel("Concentration of Br(mol/L)");
% title("third order fit");
% hold off;

% 4th order
% n=4;
% p=polyfit(dataset(:,1),dataset(:,2),n);
% fit1=p(1)*x1.^n + p(2)*x1.^(n-1) + p(3)*x1.^(n-2) + p(4)*x1.^(n-3) + p(n+1);
% fit2=p(1)*x2.^n + p(2)*x2.^(n-1) + p(3)*x2.^(n-2) + p(4)*x2.^(n-3) + p(n+1);
% figure;
% hold on;
% plot(dataset(:,1),dataset(:,2),"ko");
% plot(x1,fit1,"k-");
% xlabel("time (min)");
% ylabel("Concentration of Br(mol/L)");
% title("fourth order fit");
% hold off;

% 5th order
n=5;
p=polyfit(dataset(:,1),dataset(:,2),n);
fit1=p(1)*x1.^n + p(2)*x1.^(n-1) + p(3)*x1.^(n-2) + p(4)*x1.^(n-3) + ...
    p(5)*x1.^(n-4) + p(n+1);
fit2=p(1)*x2.^n + p(2)*x2.^(n-1) + p(3)*x2.^(n-2) + p(4)*x2.^(n-3) + ...
    p(5)*x2.^(n-4) + p(n+1);
figure;
hold on;
plot(dataset(:,1),dataset(:,2),"ko");
plot(x1,fit1,"k-");
xlabel("time (min)");
ylabel("Concentration of Br(mol/L)");
title("fifth order fit");
hold off;

%diff
df=gradient(fit2,m/h); % rate
yln=transpose(log(-df)); % ln(rate)
xln=transpose(log(fit2)); % ln([A])
figure;
plot(xln,yln,"k*");
xlabel("Ln Cbr");
ylabel("Ln rate");

%linear regression values
n = 1.539;
ln_k = 2.309;
k = exp(ln_k); % 10.0644 mol^-0.539 / (L^-0.539 * Min)

%non-linear Regrassion
beta0 = [10.0 1.5];
x = fit2;
y = df;

