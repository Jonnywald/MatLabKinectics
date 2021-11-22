%PROBLEM SET 10
%GUILHERME BERTOLA

%Problem 1

rGa = 187; %pm
rH = 120; %pm

rGa = rGa * 0.01;% A
rH = rH * 0.01;% A
dGaH = 2.60108; % A

dGaHCH3 = 2.66043; %A

DGaNH = dGaH + rGa + rH; %A

DGaCH3 = dGaHCH3 + rGa + rH; %A

dAB = (DGaCH3 + DGaNH)/2; %A
dAB = dAB * 1e-10;  %5.7008e-10 m

%molar masses
M_GaNH = 84.74;
M_GaCH3 = 84.76;

M = (M_GaCH3*M_GaNH)/(M_GaNH+M_GaCH3); %42.374999410029496


%boltzmann constant
kb = 1.38064852e-23; %m2 kg s-2 K-1

%velocity
Ur = sqrt((8*kb*T)/(M*pi));

%Steric Factor NonLinear Molecule + NonLinear Molecule -> NonLinear Complex
P = 1e-5;

%collision theory constant
Kcoll = pi*dAB^2*Ur*exp(-E/(kb*T));

% final rate constant
K = P * Kcoll;

