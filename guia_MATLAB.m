%guia rapido do matlab%

% declaracao de variaveis
a = 1;
b = a;
c = a + b;
% vetores e matrizes%
a = [1 2 3 4];
b = [1 2 3; 4 5 6; 7 8 9];
% transposta
c = b';
% indexacao de matriz
b(1,2);
%texto
t = "hello, world!";
q = "texto com ""aspas"", legal!";

%PLOTAGEM DE GRAFICOS

x = 0:pi/100:2*pi; %cira um vetor de 0 a 2pi com intervalos de um centesimo de pi
y = sin(x); %calcula o seno de cada valor de x e cria o vetor y com o resultado

plot(x,y); %plota o valor de y para cada x
xlabel("X"); %titulo do eixo X
ylabel('sin(x)');%titulo do eixo Y
title('Plot of the Sine Function'); %titulo do grafico
plot(x,y, "r--") %red dashed line

hold on %mantem o grafico aberto para se plotar por cima

y2 = cos(x); % calculo do conseno dos valores do vetor x
plot(x,y2,':') %plota os valores de y por x com a linha pontilhada
legend('sin','cos') %adiciona legenda para cada plotagem

hold off % finaliza a exibicao

%graficos 3D
[X,Y] = meshgrid(-2:.2:2);                                
Z = X .* exp(-X.^2 - Y.^2);
surf(X,Y,Z)


