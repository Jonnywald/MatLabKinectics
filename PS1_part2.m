%PS1 part 2
% Guilherme Bertola

t_range = linspace(200,400,201); %create a vector of 201 temperatures between 200K and 400K
p_range = [0.1 1.0 10.0]; %Vector of three presures 0.1 bar, 1.0 bar and 10.0 bar

%bellow is the heat capacities coheficients
a_NO2 = 22.929;
a_N2O4 = 33.054;
b_NO2 = 5.711 * 100;
b_N2O4 = 18.661 * 100;
c_NO2 = -3.519 * 100000;
c_N2O4 = -11.339 * 100000;
d_NO2 = 7.866 * 1000000000;
d_N2O4 = 0;

t0 = 298.15; %temperature of the standart state in Kelvin
R = 0.0831446261815324; %R gas constant L*bar*K-1*mol-1
dHt0 = 9.2; % heat of formation of N2O4 at 25C in KJ/Mol
dGt0 = 97.9; % gibbs free energy of N2O4 at 25C in KJ/Mol
Kat0 = exp((dGt0)/(R*t0)); %calc of the Ka at the standart state
% [NO2]
% [N2O4]
stoicMat = [2 -1]; %stoichimetric matrix of the reaction

da = stoicMat(1)*a_NO2 + stoicMat(2)*a_N2O4; %the Delta "a" of the components
db = stoicMat(1)*b_NO2 + stoicMat(2)*b_N2O4; %the Delta "b" of the components
dc = stoicMat(1)*c_NO2 + stoicMat(2)*c_N2O4; %the Delta "c" of the components
dd = stoicMat(1)*d_NO2 + stoicMat(2)*d_N2O4; %the Delta "d" of the components

lnkatkat0 = [];
kat = [];
i = 1;
for t = t_range
    lnkatkat0(i) = (da/R)*log(t/t0)+(db/(2*R))*(t-t0)+(dc/(6*R))*(t^2-t0^2)+(dd/(12*R))*(t^3-t0^3)+(1/R)*(-dHt0+da*t0+db*((t0^2)/2)+dc*((t0^3)/3)+dd*((t0^4)/4))*((1/t)-(1/t0));
    kat(i) = exp(lnkatkat0(i) + log(Kat0));
    i = i+1;
end

