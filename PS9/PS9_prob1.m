%problem set 9
%Guilherme Bertola
%problem 1

cCH40 = 2.265e18; %molecules/cm3
cO20 = 4.53e18; %molecules/cm3
cAr0 = 1.79e19; %molecules/cm3
cCO20 = 1.0; %molecules/cm3
cH20 = 1.0; %molecules/cm3
 
x0 = zeros(1,15);
x0(1) = cCH40;
x0(10) = cO20;
x0(14) = cAr0;
x0(13) = cCO20;
x0(7) = cH20;
x0(2) = 1e20;
timespan = logspace(-10,-2,200); %seconds
% timespan = linspace(0,0.01,100); %seconds
[t,x] = ode15s(@prob,timespan,x0);
figure;
loglog(t,x)
legend("CH_4","M","CH_3","H","OH","H_2O","H_2","O","CH_2O","O_2","CHO","CO","CO_2","Ar","HO_2");
xlabel("time (s)");
ylabel("Concentration");
CH4 = x(:,1);
CO2 = x(:,13);
O2 = x(:,10);
H2O = x(:,6);
figure;
loglog(t,CH4,"k-",t,CO2,"k--",t,O2,"k*",t,H2O,"ko");
legend("CH_4","CO_2","O_2","H_2O");
xlabel("time (s)");
ylabel("Concentration");

function f=prob(t,x)
    %reagents
    CH4 = x(1);
    M = x(2);
    CH3 = x(3);
    H = x(4);
    OH = x(5);
    H2O = x(6);
    H2 = x(7);
    O = x(8);
    CH2O = x(9);
    O2 = x(10);
    CHO = x(11);
    CO = x(12);
    CO2 = x(13);
    Ar = x(14);
    HO2 = x(15);
    %Forward rates const.
    k1 = 7.2e-17;
    k2 = 4.3e-11;
    k3 = 3.3e-11;
    k4 = 3.6e-12;
    k5 = 1.66e-10;
    k6 = 3.32e-14;
    k7 = 2.63e-11;
    k8 = 1.84e-10;
    k9 = 8.7e-12;
    k10 = 6.38e-16;
    k11 = 1.66e-10;
    k12 = 1.66e-10;
    k13 = 3.32e-10;
    k14 = 7e-14;
    k15 = 8.9e-13;
    k16 = 3e-11;
    k17 = 1.2e-11;
    k18 = 5.34e-12;
    k19 = 5.8e-33;
    k20 = 9.65e-32;
    k21 = 2.58e-10;
    k22 = 5.32e-33;
    k23 = 1.57e-11;
    %Reverse Rates Const.
    k1r = 3.05e-31;
    k2r = 1.74e-13;
    k3r = 1.34e-12;
    k4r = 1.1e-13;
    k5r = 8.8e-17;
    k6r = 7.2e-20;
    k7r = 1.5e-14;
    k8r = 1.42e-14;
    k9r = 7e-15;
    k10r = 5.4e-32;
    k11r = 1.9e-19;
    k12r = 2.1e-20;
    k13r = 4.6e-19;
    k14r = 1e-35;
    k15r = 3.8e-13;
    k16r = 3e-12;
    k17r = 1.1e-11;
    k18r = 2.3e-11;
    k19r = 5.14e-21;
    k20r = 8.56e-20;
    k21r = 1e-15;
    k22r = 3.6e-14;
    k23r = 2e-12;
    %forward reaction rates
    r1f = k1*CH4*M;
    r2f = k2*CH4*OH;
    r3f = k3*CH4*H;
    r4f = k4*CH4*O;
    r5f = k5*CH3*O;
    r6f = k6*CH3*O2;
    r7f = k7*CH2O*O2;
    r8f = k8*CH2O*OH;
    r9f = k9*CH2O*H;
    r10f = k10*CH2O*M;
    r11f = k11*CHO*O;
    r12f = k12*CHO*OH;
    r13f = k13*CHO*H;
    r14f = k14*CHO*M;
    r15f = k15*CO*OH;
    r16f = k16*H2*OH;
    r17f = k17*H2*O;
    r18f = k18*H*O2;
    r19f = k19*H*OH*Ar;
    r20f = k20*H*OH*H2O;
    r21f = k21*H*HO2;
    r22f = k22*H*O2*M;
    r23f = k23*OH*OH;
    %reverse reaction rates
    r1r = k1r*CH3*H*M;
    r2r = k2r*CH3*H2O;
    r3r = k3r*CH3*H2;
    r4r = k4r*CH3*OH;
    r5r = k5r*CH2O*H;
    r6r = k6r*CH2O*OH;
    r7r = k7r*CHO*OH;
    r8r = k8r*CHO*H2O;
    r9r = k9r*CHO*H2;
    r10r = k10r*CHO*H*M;
    r11r = k11r*CO*OH;
    r12r = k12r*CO*H2O;
    r13r = k13r*CO*H2;
    r14r = k14r*CO*H*M;
    r15r = k15r*CO2*H;
    r16r = k16r*H*H2O;
    r17r = k17r*H*OH;
    r18r = k18r*O*OH;
    r19r = k19r*H2O*Ar;
    r20r = k20r*H2O*H2O;
    r21r = k21r*OH*OH;
    r22r = k22r*HO2*M;
    r23r = k23r*H2O*O;
    %reaction rates
    r1 = r1f-r1r;
    r2 = r2f-r2r;
    r3 = r3f-r3r;
    r4 = r4f-r4r;
    r5 = r5f-r5r;
    r6 = r6f-r6r;
    r7 = r7f-r7r;
    r8 = r8f-r8r;
    r9 = r9f-r9r;
    r10 = r10f-r10r;
    r11 = r11f-r11r;
    r12 = r12f-r12r;
    r13 = r13f-r13r;
    r14 = r14f-r14r;
    r15 = r15f-r15r;
    r16 = r16f-r16r;
    r17 = r17f-r17r;
    r18 = r18f-r18r;
    r19 = r19f-r19r;
    r20 = r20f-r20r;
    r21 = r21f-r21r;
    r22 = r22f-r22r;
    r23 = r23f-r23r;
    
    %Mass balances
    dCH4 = -r1-r2-r3-r4;
    dM = r1+r10+r14+r22;
    dCH3 = r1+r2+r3+r4-r5-r6;
    dH = r1-r3+r5-r9+r10-r13+r14+r15+r16+r17-r18-r19-r20-r21-r22;
    dOH = -r2+r4+r6+r7-r8+r11-r12-r15-r16+r17+r18-r19-r20+r21-r23;
    dH2O = r2+r8+r12+r16+r19+r20+r23;
    dH2 = r3+r9+r13-r16-r17;
    dO = -r4-r5-r7-r11-r17+r18+r23;
    dCH2O = r5+r6-r7-r8-r9-r10;
    dO2 = -r6-r18-r22;
    dCHO = r7+r8+r9+r10-r11-r12-r13-r14;
    dCO = r11+r12+r13+r14-r15;
    dCO2 = r15;
    dAr = r19;
    dHO2 = -r21+r22;
    
    f = [dCH4;dM;dCH3;dH;dOH;dH2O;dH2;dO;dCH2O;dO2;dCHO;dCO;dCO2;dAr;dHO2];
end