%PROBLEM SET 10
%GUILHERME BERTOLA

%Problem 1

rGa = 187; %pm
rH = 120; %pm

rGa = rGa * 0.01;% A
rH = rH * 0.01;% A
dGaH = 2.47032; % A

dGaHCH3 = 2.58191; %A

DGaNH = dGaH + rGa + rH; %A

DGaCH3 = dGaHCH3 + rGa + rH; %A

dAB = (DGaCH3 + DGaNH)/2; %A

M_GaNH = 84.74;
M_GaCH3 = 84.76;

M = (M_GaCH3*M_GaNH)/(M_GaNH+M_GaCH3);
