%Problem Set 5
%Problem 2

%RAW 6.10

% A -> B

% CSTR

% Species | MolFeed Fi0| Exit (Fi) |  Exit(Fi)|       |     
%---------|------------|-----------|----------|-------|
% A       |    Fa0     |    Fa     |  Fa0 - z |       |
% B       |    0       |    Fb     |   + z    |       |
%---------|------------|-----------|----------|-------|
% Total   |Ff = Fa0    |    FT     |   Fa0    |       |

% Known Parameters

Ea = 7554.0; %K
T0 = 298.0; %K
Cp = 4.19; %Kj(Kg K)
Ca0 = 3.0; %Kmol/m3
Vr = 0.018; %M3
k0 = 4.48e6; %s-1
dHr = -2.09e5; %Kj/kmol
dens = 1e3; %Kg/M3
V0 = 6.0e-5; %m3/s

R = 8.31446261815324/1000;

% a) what T to achive Xa = 0.8

Xa = 0.8; %Conversion

Ca = Ca0 * ( 1 - Xa ); %final concentration calculated

Fa0 = Ca0*V0; %intial molar flow (Kmol/s)
Fa = Ca*V0; %final molar flow (Kmol/s)

k1 = k0*exp(-Ea/R/T0);

%mol malance
% 0 = Fa0 - Fa + r*Vr
%r = k*Ca
%energy balance
% 0 = Q - Fa0*Cp*(T-T0) - [dHr+Cp*(T-T0)]*Fa0*Xa
x0 = [300 10];
Values = fsolve(@(x)func(x,Ea,R,T0,Fa,Fa0,Vr,dHr,Cp,Xa,k1,Ca),x0);

temp = Values(1); %300.0544289272166 K
Heat = Values(2); %-30.093210972615950 KJ/s

%b) I'm cooling the reactor, because my Q (heat) is negative, meaning that
%my reaction is releasing heat, wich needs to be absorved by my
%Heating/Cooling system (cool)


function f=func(x,Ea,R,T0,Fa,Fa0,Vr,dHr,Cp,Xa,k1,Ca)
    T = x(1);
    k = k1*exp(-Ea/R * (1/T0 - 1/T));
    ra = k * Ca;
    Q = x(2);
    f(1) = Fa0 - Fa + ra*Vr;
    f(2) = Q - Fa0*Cp*(T-T0) - (dHr+Cp*(T-T0))*Fa0*Xa;

end



