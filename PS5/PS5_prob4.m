%Problem set 5
% problem 4

% Handout P 8-9 (a-e)

% 2 VinylAcetylene -> Stylene
% 2A -> B



Ca0 = 1.0; % mol/dm3
Fa0 = 5.0; % mol/s

% dHr = -231-0.012(T - 298)  %Kj/Mol

Cpa = 0.1222 %Kj/mol*K

%K = 1.48e11 * exp(-19.124/T) %dm3/mol*s

T0 = 675; %K

Ua = 5.0; %Kj/s*dm3/K


Ta = 700; %K

%a) determine the Xa in a Vr = 10 dm3 PFR with T0 = 675
Vr = 10; %dm3

function f=part_a(t,x)
    Ca = x(1);
    
    k = 1.48e11 * exp(-19.124/T);
    ra = k*Ca^2;
    dCadv = ra/V;
    f = [dCadv];
end