%problem set 8
%guilherme Bertola

%problem 4
%Non-linear Regression

Ph2 = [0.00922;0.0136;0.0197;0.028;0.0291;0.0389;0.485];
Pno = [0.00918;0.0184;0.0298;0.0378;0.0491];
rh20_Ph2 = [1.60;2.56;3.27;3.64;3.48;4.46;4.75];
rh20_Pno = [1.47;2.48;3.45;4.06;4.75];
rh20_Ph2 = rh20_Ph2./1e5;
rh20_Pno = rh20_Pno./1e5;
% val = [Ph2,Pno];
test1 = [0.05;0.05;0.05;0.05;0.05];
test = [0.05;0.05;0.05;0.05;0.05;0.05;0.05];
valh2 = [Ph2,test];
valno = [test1,Pno];
k0 = [10 10 10];
K_rate1_h2 = nlinfit(valh2,rh20_Ph2,@rate1,k0);
K_rate2_h2 = nlinfit(valh2,rh20_Ph2,@rate2,k0);
K_rate3_h2 = nlinfit(valh2,rh20_Ph2,@rate3,k0);

h2rate_1_k = K_rate1_h2(1);     %0.337362754531681
h2rate_1_kno = K_rate1_h2(2);   %0.235802546335091
h2rate_1_kh2 = K_rate1_h2(3);   %76.654556095781540

h2rate_2_k = K_rate2_h2(1);     %1.404564422445006
h2rate_2_kno = K_rate2_h2(2);   %0.440949454633404
h2rate_2_kh2 = K_rate2_h2(3);   %0.001119925813275

h2rate_3_k = K_rate3_h2(1);     %0.130775607362311
h2rate_3_kno = K_rate3_h2(2);   %0.041312237262566
h2rate_3_kh2 = K_rate3_h2(3);   %7.128203290693986

K_rate1_no = nlinfit(valno,rh20_Pno,@rate1,k0);
K_rate2_no = nlinfit(valno,rh20_Pno,@rate2,k0);
K_rate3_no = nlinfit(valno,rh20_Pno,@rate3,k0);

norate_1_k = K_rate1_no(1);     %-3.384184857084017
norate_1_kno = K_rate1_no(2);   %-0.458281230620010
norate_1_kh2 = K_rate1_no(3);   %1.436391807183324e+03

norate_2_k = K_rate2_no(1);     %0.012619645024284
norate_2_kno = K_rate2_no(2);   %17.457484573451346
norate_2_kh2 = K_rate2_no(3);   %0.008107170778276

norate_3_k = K_rate3_no(1);     %0.188143733004036
norate_3_kno = K_rate3_no(2);   %6.796863297046937
norate_3_kh2 = K_rate3_no(3);   %0.026733310354303

function r=rate1(K,val)
    k = K(1);
    Kno = K(2);
    Kh2 = K(3);
    Ph2 = val(:,1);
    Pno = val(:,2);
    r=(k.*Kno.*Pno.*Ph2)./(1+Kno.*Pno+Kh2.*Ph2);
end
function r=rate2(K,val)
    k = K(1);
    Kno = K(2);
    Kh2 = K(3);
    Ph2 = val(:,1);
    Pno = val(:,2);
    r=(k.*Kh2.*Kno.*Pno)./(1+Kno.*Pno+Kh2.*Ph2);
end
function r=rate3(K,val)
    k = K(1);
    Kno = K(2);
    Kh2 = K(3);
    Ph2 = val(:,1);
    Pno = val(:,2);
    r=(k.*Kno.*Kh2.*Pno.*Ph2)./(1+Kno.*Pno+Kh2.*Ph2).^2;
end