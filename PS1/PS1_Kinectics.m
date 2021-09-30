% Problems SET 1
% Guilherme Bertola

% 1.) RAWLINGS & EKERDT
% 2.1
%(a)
% [N2O5]
% [N2O4]
% [O2]
% [NO2]
% [NO3]
% [NO]

StoicMat = [-2 2 1 0 0 0;
            -1 0 0 1 1 0;
            1 0 0 -1 -1 0;
            0 0 1 0 -1 1;
            -1 0 0 3 0 -1;
            0 1 0 -2 0 0];
        
% (b)
rank(StoicMat);
% there are 4 linearly independent

% (c)
StoicMat2 = StoicMat;
for i=1:6
    for j=1:6
        StoicMat2([i j], : )=[];
        if rank(StoicMat2) == 4 && i~=j
            disp("remove: " + i + " & " + j);
        end
        StoicMat2 = StoicMat;
    end
end
% the common removed reactions were Eq2 and Eq3
StoicMatInd = StoicMat;
StoicMatInd([2 3], : )=[];
numrow = height(StoicMatInd);
rank(StoicMatInd);
% removing reactions 2 and 3 we are left with a rank 4, 4 reactions set

% (d) 

% no because there is no combinations of removed reactions in which both
% reactions 2 and 3 are present and the number of the rank is equal to 4 (the original rank)


% 2.15
% (a)
%[A]
%[B]
%[C]

StMat = [-1 1 0;
         0 -2 1];
Ra = -4.0;
Rb = 2.2;
Rc = 1.0;
Rm = [-4.0;2.2;1.0];
TStMat = StMat.';
%rest = lsqr(StMat,Rm)
rest = (inv(StMat*TStMat))*StMat*Rm;

r1 = rest(1);
r2 = rest(2);

% (b)


% Ra = -r1
% Rb = r1 -2r2
% Rc = r2

% (c)

RaCalc = -r1; % -4.033
RbCalc = r1 - 2*r2;% 2.1667
RcCalc = r2;% 0.933

%the values are close but no exact, since it was calculated using a
%aproximation
r = -0.25 + (0.25+0.25)*rand(3,500);
hold on
for i=1:500
    r(1,i) = (r(1,i)+1)*Rm(1,1);
    r(2,i) = (r(2,i)+1)*Rm(2,1);
    r(3,i) = (r(3,i)+1)*Rm(3,1);
    NewRm = [r(1,i);r(2,i);r(3,i)];
    rest = (inv(StMat*TStMat))*StMat*NewRm;
    r1 = rest(1,1);
    r2 = rest(2,1);
    plot(r1,r2,"k o");
    
end
xlabel("r1");
ylabel("r2");
title('Plot of estimated mesurements (with 25% deviation)');
hold off


%2.) P1-17
% (a)

y0 = [500;200]; % rabbits foxes
ts = [0 800]; % days

[t,y] = ode45(@myfunc,ts,y0);

plot(t,y)
xlabel("Days Passed")
ylabel("Population size")
legend("rabbits", "foxes")
title("Populations of rabbits and foxes by days passed");

% x1 = y(:,1);
% x2 = y(:,2);
% plot(x1,x2)
% xlabel("Rabbits")
% ylabel("Foxes")
% title("Populations of rabbits by populations of foxes");
% 
% % The curve shows that with the pre determined set of death and birth rates
% % of the rabbits and foxes, there a equilibrium between the two, rabbits
% % and foxes control their populations by the number of the other, adn with
% % a set number of rabbits you can find the highest or the lowest number of
% % foxes

%by changing the growth of rabbits to a greater value it can be noted that
%the foxes population quickly catches up to the rabbit population.
%by doublig the value it can be noted that the foxes population have the
%same behavior as the rabbit, up to that value the population of foxes then
%quickly becomes unstable.
%and as the death rate of the rabbit is decresed it is possible to see how
%the population of of foxes quickly surpasses the rabbits, the death rates
%should be thinkerd with caution, since it has a bigger impact than the
%growth rate of both rabbits and foxes







