% PROBLEM SET 3
% Guilherme Bertola
% Problem 3

% A -> P

% Species  | Feed | Exit (Ni)   |       yi         |     
%--------------------------------------------------|
% A        | Ca0  |    Ca3      | (Ca3)/(Ca3 + Cp) |
% P        |  0   |    Cp       | (Cp)/(Ca3 + Cp)  |
%----------|------|-------------|------------------|
% Total    | Ca0  |   Ca3 + Cp  |        1         |

Ca0 = 1.0; %Mol / L
Ca1_0= 0.0;%Mol / L
Ca2_0= 0.0;%Mol / L
Ca3_0= 0.0;%Mol / L

k1 = 1.0; %min-1
k2 = 2.0; %min-1
k3 = 5.0; %min-1

V0 = 1000.0; %L/min
Vr = 100.0; %L/min

V1 = 1000.0; % L
V2 = 1500.0; % L
V3 = 100.0; % L

tau1 = V1/V0; %min
t_span = linspace(0,5*tau1,200); %time span of the reaction o to 5 times tau1

% the call to the ode45
[y,C] = ode45(@(t,x)prob3(t,x,V0,Vr,V1,V2,V3,k1,k2,k3,Ca0),t_span,[Ca1_0 Ca2_0 Ca3_0]);

%concentrations at each reactor
Ca1 = C(:,1);
Ca2 = C(:,2);
Ca3 = C(:,3);

%Mol flow after each reactor
Fa1 = Ca1.*(V0+Vr);
Fa2 = Ca2.*(V0+Vr);
Fa3 = Ca3.*Vr;

%Mol flow before each reactor
Fa1_0 = V0.*Ca0 + Vr.*Ca2; 
Fa2_0 = Ca1.*(V0+Vr);
Fa3_0 = Vr.*Ca2;

%cumulative conversion after each reactor
Xa1 = (Fa1_0-Fa1)./Fa1_0;
Xa2 = (Fa1_0-Fa2)./Fa1_0;
Xa3 = (Fa1_0-Fa3)./Fa1_0;

hold on;
plot(t_span,Xa1,"- k");
plot(t_span,Xa2,"-- k");
plot(t_span,Xa3,": k");
xlabel("time (min)");
ylabel("Conversion");
title("Cumulative conversion of A in each reactor");
legend("Xa1","Xa2","Xa3");
hold off;

Ca0 = 0.0;
% the call to the ode45
[y,C] = ode45(@(t,x)prob3(t,x,V0,Vr,V1,V2,V3,k1,k2,k3,Ca0),t_span,[Ca1_0 Ca2_0 Ca3_0]);

%concentrations at each reactor
Ca1 = C(:,1);
Ca2 = C(:,2);
Ca3 = C(:,3);

%Mol flow after each reactor
Fa1 = Ca1.*(V0+Vr);
Fa2 = Ca2.*(V0+Vr);
Fa3 = Ca3.*Vr;

%Mol flow before each reactor
Fa1_0 = V0.*Ca0 + Vr.*Ca2; 
Fa2_0 = Ca1.*(V0+Vr);
Fa3_0 = Vr.*Ca2;

%cumulative conversion after each reactor
Xa1 = (Fa1_0-Fa1)./Fa1_0;
Xa2 = (Fa1_0-Fa2)./Fa1_0;
Xa3 = (Fa1_0-Fa3)./Fa1_0;


% hold on;
% plot(t_span,Xa1,"- k");
% plot(t_span,Xa2,"-- k");
% plot(t_span,Xa3,": k");
% xlabel("time (min)");
% ylabel("Conversion");
% title("Cumulative conversion of A in each reactor");
% legend("Xa1","Xa2","Xa3");
% hold off;

function f = prob3(t,x,V0,Vr,V1,V2,V3,k1,k2,k3,Ca0)
    dCa1dt = (V0/V1)*Ca0 + (Vr/V1)*x(2) - ((V0+Vr)/V1)*x(1) - k1*x(1);
    dCa2dt = ((V0+Vr)/V2)*x(1) - ((V0+Vr)/V2)*x(3) - k2*x(2);
    dCa3dt = (V0/V3)*x(2) - (V0/V3)*x(3) - k3*x(3);
    f = [dCa1dt;dCa2dt;dCa3dt];
end