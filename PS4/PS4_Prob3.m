%Problem set 4
% Problem 3

% M-Xylene -> benzene + methane (1)
% M-xylene -> p-xylene (2)

% mx = M-xylene
% ben = Benzene
% met = Methane
% px = P-xylene

% Species | MolFeed Fi0| Exit (Fi) |    yi    |     
%---------|------------|-----------|----------|
% mx      |    Fmx     |    Fmx    |  Fmx/FT  |
% ben     |     0      |    Fben   |  Fben/FT |
% met     |     0      |    Fmet   |  Fmet/FT |
% px      |     0      |    Fpx    |  Fpx/FT  |
%---------|------------|-----------|----------|
% Total   |    Fmx     |    FT     |     1    |

%question A
Xmx_obj = 0.9;
k1 = 0.22; %s-1
k2 = 0.71; %s-1
T = 673; % C
ymx0 = 0.75; %composition at the feed (25% inerts)
V0 = 2000; %l/min volumetric feed
Ct0 = 0.05; %mol/l total concentration of the feed

%calculating Cmx0
Cmx0 = Ct0 * ymx0;

%converting the volumetric feed to liter per seconds
V0 = V0 / 60;

%creating a vector of volumes
v_span = linspace(0,100,500);

[v,c] = ode45(@(v,c)func(v,c,V0,k1,k2),v_span,[Cmx0 0 0 0]);

%assigning the concentraions
Cmx = c(:,1);
Cben = c(:,2);
Cmet = c(:,3);
Cpx = c(:,4);

%calculating the convertions
Xmx = (Cmx0 - Cmx)./Cmx0 ;

%finding the volume that Xmx >= 0.9
i = 1;
for x=Xmx.'
    if (x >= 0.9)
        break
    end
    i = i+1;
end
vol = v(i);

%calculating the space-time
tau = vol / V0; % 2.476953907815631 Seconds

vflow = linspace(0,V0,100);

taus = vol ./ vflow;
range = linspace(1,100,100);
yield = zeros(1,100);
sel = zeros(1,100);
for i=range
    sel(i) = (Cpx(i)*vflow(i))/((Cmx(i)+Cben(i)+Cmet(i))*vflow(i));
    yield(i) = (Cpx(i)*vflow(i))/(((Cmx0)*vflow(i))-((Cmx(i)+Cben(i)+Cmet(i))*vflow(i)));
end


hold on
plot(taus,sel);
plot(taus,yield);
set(gca, 'XDir','reverse');
xlabel("Space-time 'tau' (s)");
ylabel("Yield and Selectivity");
legend("Selectivity","Yield");
title("Yield and Selectivity by Space-time");
hold off

%question B (b)

tau2 = 0.5; %s
E1 = 20000; %cal/mol
E2 = 10000; %cal/mol
R = 1.98720425864083; %cal/K*mol

Ts = 946.15; % K

%calculating the k0 of each reaction
k1_0 = k1/exp(-E1/(R*Ts));
k2_0 = k2/exp(-E2/(R*Ts));

%calculating the volume of the CSTR
Vcstr = tau2 * V0; 

%inital concentrations
Cmx0 = Ct0 * ymx0;
Cpx0 = 0;
Cben0 = 0;
Cmet0 = 0;
Cit = Ct0 * (1 - ymx0);

%range of temperatures in kelvin
t_range = linspace(300,3000,1000); %K

Cpx_2 = zeros(1,1000);
Cmx_2 = zeros(1,1000);
i = 1;

for t=t_range
    
    k1T = k1_0*exp(-E1/(R*t));
    k2T = k2_0*exp(-E2/(R*t));
    Cmx_2(i) = Cmx0 / (1 + k1T*tau2 +k2T*tau2);
    Cpx_2(i) = (k2T*Cmx_2(i))*tau2;
    i = i + 1;
end

maxC = max(Cpx_2); %0.013072912471449 mol/L
val = find(Cpx_2==maxC,1);

BestTemp = t_range(val); %1194.59 K


function f = func(v,c,V0,k1,k2)
dmxdv = (-k1*c(1) - k2*c(1))/V0;
dbendv = (k1*c(1))/V0;
dmetdv = (k1*c(1))/V0;
dpxdv = (k2*c(1))/V0;
f = [dmxdv;dbendv;dmetdv;dpxdv];
end