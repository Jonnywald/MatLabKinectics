%PS1 part 2
% Guilherme Bertola

t_range = linspace(200,400,500); %create a vector of 201 temperatures between 200K and 400K
p_range = [0.1 1.0 10.0]; %Vector of three presures 0.1 bar, 1.0 bar and 10.0 bar
V = 100; % 100 L of volume
%bellow is the heat capacities coheficients
a_NO2 = 22.929;
a_N2O4 = 33.054;
b_NO2 = 5.711 / 100;
b_N2O4 = 18.661 / 100;
c_NO2 = -3.519 / 100000;
c_N2O4 = -11.339 / 100000;
d_NO2 = 7.866 / 1000000000;
d_N2O4 = 0;

t0 = 298.15; %temperature of the standart state in Kelvin
R = 0.0831446261815324; %R gas constant L*bar*K-1*mol-1
dHt0 = 33.2; % heat of formation of NO2 at 25C in KJ/Mol
dGt0 = 51.3; % gibbs free energy of NO2 at 25C in KJ/Mol
Kat0 = exp((-dGt0)/(R*t0)); %calc of the Ka at the standart state
% [NO2]
% [N2O4]
stoicMat = [2 -1]; %stoichimetric matrix of the reaction

da = stoicMat(1)*a_NO2 + stoicMat(2)*a_N2O4; %the Delta "a" of the components
db = stoicMat(1)*b_NO2 + stoicMat(2)*b_N2O4; %the Delta "b" of the components
dc = stoicMat(1)*c_NO2 + stoicMat(2)*c_N2O4; %the Delta "c" of the components
dd = stoicMat(1)*d_NO2 + stoicMat(2)*d_N2O4; %the Delta "d" of the components

%
%  Ka(t) = (aNO2)2
%          -------
%          (aN2O4)
%

% NiNO2 = 2*Zeta
% NiN2O4 = Ni0 - Zeta
% Total = Ni0 + zeta

%yNO2 = (2*Zeta) / (Ni0 + zeta)
%yN2O4 = (Ni0 - Zeta) / (Ni0 + Zeta)

% Pt = (N*R*T)/V
%
% Pt = ((Ni0+zeta)*R*T)/V

% aNO2 = (2*Zeta)   (1+Zeta)RT
%        -------- * -----------
%        (1+Zeta)       V

% aN2O4 =(1-Zeta)   (1+Zeta)RT
%        -------- * -----------
%        (1+Zeta)       V

lnkatkat0 = []; % matrix for ln of Kat/kat0
CHR_lnkatkat0 = []; % constant heat of reaction ln of Kat/kat0
kat = [];% matirx for Kat
CHR_Kat = []; %constant heat of reaction Kat
i = 1;
j = 1;
zeta = []; % matrix for the extent of reaction
CHR_zeta = []; % constant heat of reaction extent of reaction
x0 = 1.0001; % inital guess for fsolve
global pt;
global k;
for p = p_range
    for t = t_range
        lnkatkat0(i,j) = (da/R)*log(t/t0)+(db/(2*R))*(t-t0)+(dc/(6*R))*(((t).^2)-((t0).^2))+(dd/(12*R))*(((t).^3)-((t0).^3))+(1/R)*((-dHt0)+(da*t0)+(db*((t0.^2)/2))+(dc*((t0.^3)/3))+(dd*((t0.^4)/4)))*((1/t)-(1/t0));
        CHR_lnkatkat0(i,j) = (-dHt0/R)*((1/t)-(1/t0));
        kat(i,j) = exp(lnkatkat0(i,j) + log(Kat0));
        CHR_Kat(i,j) = exp(CHR_lnkatkat0(i,j) + log(Kat0));
        pt = p;
        k = kat(i,j);
        zeta(i,j) = fsolve(@myfunc,x0);
        k = CHR_Kat(i,j);
        CHR_zeta(i,j) = fsolve(@myfunc,x0);
        i = i+1;
    end
    i = 1;
    j = j + 1;
end

hold on
plot(t_range,zeta)
xlabel("temperature (K)")
ylabel("Extent of reaction")
legend("0.1 bar","1.0 bar","10.0 bar")
title("Extent of reaction over temperature (varying pressure)")

plot(t_range,CHR_zeta)
xlabel("temperature (K)")
ylabel("Extent of reaction")
legend("0.1 bar","1.0 bar","10.0 bar")
title("Extent over temp (varying pressure & constant heat of reaction)")
hold off
function u = myfunc(z)
    global pt;
    global k;
    u = (((((2*z)/(1+z))*pt).^2)/(((1-z)/(1+z))*pt)) - k;
end






