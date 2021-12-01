%PROBLEM SET 10
%GUILHERME BERTOLA

%Problem 2

% A = CH
% B = N2
% T = HCN2

T = 298;
h = 6.62607004e-34;
kb = 1.38064852e-23;
c = 299792458;

MA = 13.019;
MB = 28.013;
MT = 41.032;

IA = 1.935/10000;
IB = 13.998/10000;
IT = 73.2/10000;

E0T = 221;

vA = 2733;
vB = 2330;
vT = [3130; 2102; 1252; 1170; 564; 401];

qeA = 1;
qeB = 1;
qeT = 1;

qTvA = ((2*pi*MA*kb*T)/(h^2))^(3/2);
qTvB = ((2*pi*MB*kb*T)/(h^2))^(3/2);
qTvT = ((2*pi*MT*kb*T)/(h^2))^(3/2);

qvA = 1/(1-exp((-h*c*vA)/(kb*T)));
qvB = 1/(1-exp((-h*c*vB)/(kb*T)));
qvT = 1;
for i=1:6
    qvT = qvT * (1/(1-exp((-h*c*vT(i))/(kb*T))));
end

qrA = (8*pi^2*kb*IA)/(h^2);
qrB = (8*pi^2*kb*IB)/(h^2);
qrT = (8*pi^2*kb*IT)/(h^2);

QA = qeA * qTvA * qvA * qrA;
QB = qeB * qTvB * qvB * qrB;
QT = qeT * qTvT * qvT * qrT;

Ktst = ((kb*T)/h)*(QT/(QA*QB))*exp(-E0T/(kb*T));
