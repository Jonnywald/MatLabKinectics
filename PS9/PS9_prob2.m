%problem set 9
%Guilherme Bertola
%problem 2

%CHO + O -> CO + OH

M_O2 = 15.999;
M_CHO = 29.01804;

M = (M_O2*M_CHO)/(M_O2+M_CHO); %10.312975308016698

k = 1.38064852e-23;

T = 2000;

ur = sqrt((8*k*T)/(M*pi)); %8.257234754986269e-11

d_O2 = 2*152e-12;
d_CHO = 2*235.66e-12; 
dab = (d_O2 + d_CHO)/2; %3.876600000000000e-10

e = fsolve(@(x)func(x,dab,ur,k,T),1e-18); %-1.184468525547849e-18

K_2000 = pi*dab^2*ur*exp((-e)/(k*T)); %1.660000000000004e-10

function f=func(x,dab,ur,k,T)
    f = -log((1.66e-10)/(pi .* dab.^2 .* ur))+ (-x)./(k.*T);
end