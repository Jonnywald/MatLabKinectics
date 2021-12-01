%PROBLEM SET 10
%GUILHERME BERTOLA

%Problem 3

% A = H2CO
% T = TS

NA = 6.02214076e23; % mol-1
h = 6.62607004e-34; %m2 kg / s
kb = 1.38064852e-23; %m2 kg s-2 K-1
c = 299792458; %m / s
c = c * 100; %cm / s

MA = 30.031; %g/mol
MT = 30.031; %g/mol

IAz = 6.33171; %amu * Bohr^2
IAy = 46.71640; %amu * Bohr^2
IAx = 53.04811; %amu * Bohr^2

ITz = 7.17351; %amu * Bohr^2
ITy = 51.09245; %amu * Bohr^2
ITx = 58.26594; %amu * Bohr^2

IAz = IAz * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2
IAy = IAy * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2
IAx = IAx * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2

ITz = ITz * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2
ITy = ITy * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2
ITx = ITx * 1.66054e-27 * (5.29177e-11)^2; %Kg * m-2

vA = [1198.5963; 1279.6691; 1563.0192; 1849.4365; 2916.4495; 2967.4670]; %cm-1
vT = [782.5773; 918.4587; 1285.9031; 1936.0523; 3212.5622]; %cm-1

E_A = -114.473648; %Hartree

E_T = -114.339439; %Hartree

Ea = E_T - E_A;

Ea = (Ea * 627.5095) / NA; % Kcal
Ea = Ea * 4184; %J

Temp_span = 300:100:1000;

qvA = 1;
qvT = 1;
rates = zeros(8,1);
j = 1;
for T=Temp_span
    qTvA = ((2*pi*MA*kb*T)/(h^2))^(3/2);
    qTvT = ((2*pi*MT*kb*T)/(h^2))^(3/2);
    qrA = pi^0.5 * ((8*(pi^2)*kb*T)/(h^2))^1.5 * IAx^0.5 * IAy^0.5 * IAz^0.5;
    qrT = pi^0.5 * ((8*(pi^2)*kb*T)/(h^2))^1.5 * ITx^0.5 * ITy^0.5 * ITz^0.5;
    for i=1:6
       qvA = qvA * (1/(1-exp((-h*c*vA(i))/(kb*T)))); 
    end
    for i=1:5
       qvT = qvT * (1/(1-exp((-h*c*vT(i))/(kb*T)))); 
    end
    Ktst = (kb*T)/(h) * (qTvT/qTvA) * (qrT/qrA) * (qvT/qvA) * exp((-Ea)/(kb*T));
    rates(j) = Ktst;
    j = j +1 ;
end

loglog(Temp_span,rates,"--k");
xlabel("Temperature (K)");
ylabel("Rate Constant (1/s)");
title("Rate constant vs. Temperature");


