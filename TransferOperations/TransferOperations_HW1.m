%Transfer Operations 1 - Homework Set 1

%Guilherme Bertola

%problem 1 - stress tensor

area = 1; %mm2
area = area / 1e+6; %m2

f1 = 2;  %N
f2 = -5; %N
f3 = -4; %N

p1 = f1/area; %Pa
p2 = f2/area; %Pa
p3 = f3/area; %Pa

p1 = p1 / 1e+6; %MPa
p2 = p2 / 1e+6; %MPa
p3 = p3 / 1e+6; %MPa

%a)
%state of stress at P
sigma = [p1 0 0;
         0 0 p3;
         0 p2 0];

% | 2  0  0 |
% | 0  0 -4 |
% | 0 -5  0 |

%b)
%magnitute of the force
n = [(1/sqrt(2)) 0 (1/sqrt(2))];
s = n * sigma;
s = s.*s;
s = sqrt(sum(s)); %Mpa
s = s * 1e+6; %Pa
f = s * area; %N
%F = 3.807886552931953 N

%c)
% the component of this force in N normal to the surface
s = n * sigma;
sn = s .* n;
sn = sum(sn); %MPa
sn = sn * 1e+6; %Pa
fn = sn * area; %N
% Fn = 1N - normal to the surface

% Problem 2 - torsional rheometers

%(a) shear rate expressions for each fixture

% Shear rate = gamma-dot
% radial position = r
% gap = h
% rotation speed = omega
% cone angle = alfa

% Parallel-plate fixture

%              r*omega
% gamma-dot = ---------
%                 h

% Cone-plate fixture

%              omega*r        omega*r
% gamma-dot = ---------- = -------------
%              h(r,alfa)    r*tan(alfa)

% the parallel-plate fixture has non-homogeneous shear rate field

%(b)

% the torque expression a function of r:

% M = integral(2*pi*tau*r^2,0,R);

%since the shear rate is constant in a cone fixture, the shear stress is
%also constant

%making the torque equation be (after integral):

% M = (2*pi*tau)/3 * R^3;

% isolating tau:

% tau = 3*M / 2*pi*R^3;

%(c)

% if the gap between the inner and outer cylinders of a "cup and bob"
% rheometer is too large, the flow field will no longer be "laminar" and
% vortices will appear (Taylor`s vortices), make the torque mesuared no
% compatible with the viscosity of the fluid

%Problem 3

Eta0 = 3; %Pa * s
EtaInf = 1; %Pa * s
n = 0.6;
k1 = 0.02; %s^-0.4
k2 = 0.016; %s^-1
t1 = 5.0; %s
gamma = 10.0; % s-1

%(a)

a = -(k1*(gamma^n)+k2);

Eta1 = ((-k1*(gamma^n)*exp(a*t1)-k2)/a)*(Eta0 - EtaInf)+EtaInf; %2.367094453263016 Pa * S
Tau1 = Eta1*gamma; % 23.670944532630160 Pa

%(b)
% "Second" Phase (0 - 10 s)
t2 = 10.0; % s
Eta2 = ((-k1*(gamma^n)*exp(a*t2)-k2)/a)*(Eta0 - EtaInf)+EtaInf; %1.974721018824596 Pa * S

% "Third" Phase (10 - 15 s)
t3 = 15.0; % s
m2m0 = (Eta2 - EtaInf)/(Eta0 - EtaInf);
gamma = 0;
Eta3 = (((m2m0-1)*exp(-k2*(t3-t2)))+1)*(Eta0 - EtaInf)+EtaInf; % 2.053548212870349 Pa * S
Tau3 = Eta3 * gamma; % 0 Pa

