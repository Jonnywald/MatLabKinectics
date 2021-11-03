%problem set 8
%guilherme Bertola

%problem 3
%Initial Rate method of Data analysis

dataset = [6 20 0.420; 8 20 0.647; 10 20 0.895; 12 20 1.188;...
    16 20 1.811; 10 10 0.639; 10 20 0.895; 10 40 1.265; 10 60 1.550;...
    10 100 2.021];
Pa = dataset(:,1);
Pa_const_Pb = Pa(1:5);
Pb = dataset(:,2);
Pb_const_Pa = Pb(6:end);
minus_ra = dataset(:,3);
m_ra_pa = minus_ra(1:5);
m_ra_pb = minus_ra(6:end);

beta0 = [10 1];
alpha0 = [10 1];
alphaPa = nlinfit(Pa_const_Pb,m_ra_pa,@funcPa,alpha0);
betaPb = nlinfit(Pb_const_Pa,m_ra_pb,@funcPb,beta0);
new_m_ra_pa = funcPa(alphaPa,Pa_const_Pb);
new_m_ra_pb = funcPb(betaPb,Pb_const_Pa);

figure;
hold on;
plot(Pa_const_Pb,m_ra_pa,"ko");
plot(Pa_const_Pb,new_m_ra_pa,"k-");
xlabel("Pa (torr)");
ylabel("-ra (torr/s)");
title("Pressure Pa by rate (constant Pb)");
hold off;
figure;
hold on;
plot(Pb_const_Pa,m_ra_pb,"k*");
plot(Pb_const_Pa,new_m_ra_pb,"k--");
xlabel("Pb (torr)");
ylabel("-ra (torr/s)");
title("Pressure Pb by rate (constant Pa)");
hold off;
% xlabel

alpha = alphaPa(2); % 1.488603026866449
beta = betaPb(2); % 0.504081880837406
k = nlinfit(Pa_const_Pb,m_ra_pa,@rate,10); % 0.006457432874974 torr^-0.992684907703855 * s^-1

function m_ra=funcPa(alpha,P)
    m_ra = alpha(1).*(P.^alpha(2));
end
function m_ra=funcPb(beta,P)
    m_ra = beta(1).*(P.^beta(2));
end
function m_ra=rate(val,P)
    alpha = 1.488603026866449;
    beta = 0.504081880837406;
    m_ra = val.*(P.^alpha).*(20.^beta);
end