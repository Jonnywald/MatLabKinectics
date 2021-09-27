%PRoblem set 4
% problem 2

% A + B -> C    (1)
% 2A -> D       (2)

% Species | MolFeed Fi0| Exit (Fi) |   yi    |     
%---------|------------|-----------|---------|
% A       |    FA0     |    FA     |  FA/FT  |
% B       |    FB0     |    FB     |  FB/FT  |
% C       |     0      |    FC     |  FC/FT  |
% D       |     0      |    FD     |  FD/FT  |
%---------|------------|-----------|---------|
% Total   | FA0 + FB0  |    FT     |    1    |

% r1 = k1*Ca*Cb
% r2 = k2*Ca^2

% question (a)

% higher yield is configuration 1, as this config all of A is in direct
% "contact" with reagent B meaning you are producing as much of C as you
% can, while in the config 2, a "fresh" stream of A is feeded to the second
% reactor, meaning it has "less b" to react with, making it favorable to
% generate D.

% since conversion, does not "care" to what the reagent has been converted
% it seems that the configuration 2, can achive a higher convertion of A,
% because of the favorability of the generation of D (wich consumes 2 A).

% question (b)

V1 = 1000;      %L
V2 = 1000;      %L
Vt0 = 1000;     %L/h
k1 = 1;         %L/Mol*h
k2 = 1;         %L/Mol*h
Ca0 = 2;        %Mol/L
Cb0 = 2;        %Mol/L
Cc0 = 0;        %Mol/L
Cd0 = 0;        %Mol/L

%CSTR at a Standard state

%Configuration 1

%reactor 1
% 0 = Ca0 + (-r1 -r2) * V1/Vt0  - Ca1   (1)
% 0 = Cb0 + (-r1) * V1/Vt0  - Cb1       (2)
% 0 = Cc0 + (r1) * V1/Vt0  - Cc1        (3)
% 0 = Cd0 + (r2) * V1/Vt0  - Cd1        (4)
%reactor 2
% 0 = Ca1 + (-r1 -r2) * V2/Vt0  - Ca2   (5)
% 0 = Cb1 + (-r1) * V2/Vt0  - Cb2       (6)
% 0 = Cc1 + (r1) * V2/Vt0  - Cc2        (7)
% 0 = Cd1 + (r2) * V2/Vt0  - Cd2        (8)

% configuration 2
%reactor 1
% 0 = Ca0 + (-r1 -r2) * V1/Vt0  - Ca1   (1)
% 0 = Cb0 + (-r1) * V1/Vt0  - Cb1       (2)
% 0 = Cc0 + (r1) * V1/Vt0  - Cc1        (3)
% 0 = Cd0 + (r2) * V1/Vt0  - Cd1        (4)
%reactor 2
% 0 = Ca1 + Ca0 + (-r1 -r2) * V2/Vt0  - Ca2   (5)
% 0 = Cb1 + (-r1) * V2/Vt0  - Cb2       (6)
% 0 = Cc1 + (r1) * V2/Vt0  - Cc2        (7)
% 0 = Cd1 + (r2) * V2/Vt0  - Cd2        (8)

% question (c)

% configuration 1
x0 = [1 1 1 1 0.5 0.5 1.5 1.5];
result = fsolve(@(x)config_1(x,V1,V2,Vt0,Ca0,Cb0,Cc0,Cd0,k1,k2),x0);

%Assining results
Ca1 = result(1); %0.695620769559862 Mol/L
Cb1 = result(2); %1.179509024602917 Mol/L
Cc1 = result(3); %0.820490975397083 Mol/L
Cd1 = result(4); %0.483888255043055 Mol/L
Ca2 = result(5); %0.314503336946571 Mol/L
Cb2 = result(6); %0.897303940933064 Mol/L
Cc2 = result(7); %1.102696059066949 Mol/L
Cd2 = result(8); %0.582800603986482 Mol/L

%calculating the conversion
Xa_config_1 = (Ca0-Ca2)/Ca0; %0.842748331526714
%Calculating the Yield
Yc_config_1 = Cc2/(Ca0-Ca2); %0.654226189370981

% configuration 2
result2 = fsolve(@(x)config_2(x,V1,V2,Vt0,Ca0,Cb0,Cc0,Cd0,k1,k2),x0);

%Assining results
Ca1_2 = result2(1); %0.695620769559862 Mol/L
Cb1_2 = result2(2); %1.179509024602917 Mol/L
Cc1_2 = result2(3); %0.820490975397083 Mol/L
Cd1_2 = result2(4); %0.483888255043055 Mol/L
Ca2_2 = result2(5); %1.031867081978486 Mol/L
Cb2_2 = result2(6); %0.580505011870948 Mol/L
Cc2_2 = result2(7); %1.419494988128923 Mol/L
Cd2_2 = result2(8); %1.548637929892453 Mol/L

%calculating the conversion
Xa_config_2 = (Ca0-Ca2_2)/Ca0; %0.484066459010757
%Calculating the Yield
Yc_config_2 = Cc2_2/(Ca0-Ca2_2); %1.466219112794778

function u = config_1(x,V1,V2,Vt0,Ca0,Cb0,Cc0,Cd0,k1,k2)

    Ca1 = x(1);
    Cb1 = x(2);
    Cc1 = x(3);
    Cd1 = x(4);
    Ca2 = x(5);
    Cb2 = x(6);
    Cc2 = x(7);
    Cd2 = x(8);
    
    
    u(1) = Ca0 + (-(k1*Ca1*Cb1)-(k2*Ca1^2)) * V1/Vt0  - Ca1;
    u(2) = Cb0 + (-(k1*Ca1*Cb1)) * V1/Vt0  - Cb1;
    u(3) = Cc0 + (k1*Ca1*Cb1) * V1/Vt0  - Cc1;
    u(4) = Cd0 + (k2*Ca1^2) * V1/Vt0  - Cd1;
    u(5) = Ca1 + (-(k1*Ca2*Cb2)-(k2*Ca2^2)) * V2/Vt0  - Ca2;
    u(6) = Cb1 + (-(k1*Ca2*Cb2)) * V2/Vt0  - Cb2;
    u(7) = Cc1 + (k1*Ca2*Cb2) * V2/Vt0  - Cc2;
    u(8) = Cd1 + (k2*Ca2^2) * V2/Vt0  - Cd2;

end
function u = config_2(x,V1,V2,Vt0,Ca0,Cb0,Cc0,Cd0,k1,k2)

    Ca1 = x(1);
    Cb1 = x(2);
    Cc1 = x(3);
    Cd1 = x(4);
    Ca2 = x(5);
    Cb2 = x(6);
    Cc2 = x(7);
    Cd2 = x(8);
    
    
    u(1) = Ca0 + (-(k1*Ca1*Cb1)-(k2*Ca1^2)) * V1/Vt0  - Ca1;
    u(2) = Cb0 + (-(k1*Ca1*Cb1)) * V1/Vt0  - Cb1;
    u(3) = Cc0 + (k1*Ca1*Cb1) * V1/Vt0  - Cc1;
    u(4) = Cd0 + (k2*Ca1^2) * V1/Vt0  - Cd1;
    u(5) = Ca1 + Ca0 + (-(k1*Ca2*Cb2)-(k2*Ca2^2)) * V2/Vt0  - Ca2;
    u(6) = Cb1 + (-(k1*Ca2*Cb2)) * V2/Vt0  - Cb2;
    u(7) = Cc1 + (k1*Ca2*Cb2) * V2/Vt0  - Cc2;
    u(8) = Cd1 + (k2*Ca2^2) * V2/Vt0  - Cd2;

end

