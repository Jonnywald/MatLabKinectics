%problem set 5
% problem 3

% RAW 6.14

% A <-> B
% PFR adiabatic

% Species | MolFeed Fi0| Exit (Fi) |  Exit    |     Ci      |  Ci(Xa)  |
%---------|------------|-----------|----------|-------------|----------|
% A       |    Fa0     |    Fa     |  Fa0 - z | (Fa0 - z)/V | Ca0(1-Xa)|
% B       |    Fb0     |    Fb     |  Fb0 + z | (Fb0 - z)/V | Cb0-Ca0Xa|
%---------|------------|-----------|----------|-------------|----------|
% Total   |Ff = Fa0+Fb0|    FT     |    Ff    |   Ff/V      |          |

% a) Fa = Ca(v)*V(v) , Fb = Cb(v)*V(v)
% based in the stoichiometry of the reaction we have 1:1 proportion, as
% much of Ca is consumed, as much of Cb is made, so Cb = Ca0 - Ca and the
% reverse is also true

%b) dFadv = ra  ->  dFadv = -kCa + krCb
% dFadv = -kCa + kr(Ca0 - Ca)
% dCadv = (-kCa + kr(Ca0 - Ca))/V

% V = V0 * Ft/Ff * T/T0;

%c) 
% V*ro*Cp*dTdv = -dHr*r+2/R * Q(Ta-T)
