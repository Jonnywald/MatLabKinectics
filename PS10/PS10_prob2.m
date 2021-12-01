%PROBLEM SET 10
%GUILHERME BERTOLA

%Problem 2

% A = CH
% B = N2
% T = HCN2

T = 298; % K
h = 6.62607004e-34; %m2 kg / s
kb = 1.38064852e-23; %m2 kg s-2 K-1
c = 299792458 * 100; %cm / s

MA = 13.019 * 1/1000; %kg/mol
MB = 28.013 * 1/1000; %kg/mol
MT = 41.032 * 1/1000; %kg/mol

MA = MA/NA; %Kg
MB = MB/NA; %Kg
MT = MT/NA; %Kg

IA = (1.935/10000 * 1/1000)*1e-40; %kg m2
IB = (13.998/10000 * 1/1000)*1e-40; %kg m2
IT = (73.2/10000 * 1/1000)*1e-40; %kg m2

E0T = 221; %kJ/mol
NA = 6.02214076e23; % mol-1
E0T = E0T/NA; %J

vA = 2733; %cm-1
vB = 2330; %cm-1
vT = [3130; 2102; 1252; 1170; 564; 401]; %cm-1

qeA = 1; %unitless
qeB = 1; %unitless
qeT = 1; %unitless

qTvA = ((2*pi*MA*kb*T)/(h^2))^1.5; %m-3
qTvB = ((2*pi*MB*kb*T)/(h^2))^1.5; %m-3
qTvT = ((2*pi*MT*kb*T)/(h^2))^1.5; %m-3

qvA = 1/(1-exp((-h*c*vA)/(kb*T))); %unitless
qvB = 1/(1-exp((-h*c*vB)/(kb*T))); %unitless
qvT = 1;
for i=1:6
    qvT = qvT * (1/(1-exp((-h*c*vT(i))/(kb*T)))); %unitless
end

qrA = (8*pi^2*kb*IA*T)/(h^2);%unitless
qrB = (4*pi^2*kb*IB*T)/(h^2);%unitless %homonuclear
qrT = (8*pi^2*kb*IT*T)/(h^2);%unitless

QA = qeA * qTvA * qvA * qrA; %m-3
QB = qeB * qTvB * qvB * qrB; %m-3
QT = qeT * qTvT * qvT * qrT; %m-3

Ktst = ((kb*T)/h)*(QT/(QA*QB))*exp(-E0T/(kb*T)); % m3/s
Ktst = Ktst * NA * 1e3; %1.227e+08 L/(mol*s)
