clear all;
clear;
clc;

% Depósito de combustible en perfil alar. Apartado 1

% Perfil %
% Fuente y borde de ataque
N = 5; % Densidad de mallado
xc = [0; 0; 0.2; 0.2];
yc = [0; .2; 0.2; 0];
[x y tri F] = mesh_block(xc, yc, 2*N, 2*N);

xc = [0; 0; 0.2; 0.2];
yc = [0; -0.1; -0.1; 0];
dc = [0; -1; 1; 0];
[y1 x1 tri1 F1] = mesh_block_hermite(xc, yc, dc, 2*N, 2*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [1.2; 1.2; 1.6; 1.6];
yc = [0; .2; 0.2; 0];
[x1 y1 tri1 F1] = mesh_block(xc, yc, 2*N, 4*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [1.6; 1.6; 2; 2];
yc = [0; 0.2; 0.1; 0];
dc = [0.25; -0.25; -0.25; 0];
[x1 y1 tri1 F1] = mesh_block_hermite(xc, yc, dc, 2*N, 4*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [2; 2; 2.4];
yc = [0; 0.1; 0];
[x1 y1 tri1 F1] = mesh_triangle(xc, yc, 2*N, 4*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

% Succion
xc = [0; 0; -0.1];
yc = [0.2; 0.3; 0.2];
[x1 y1 tri1 F1] = mesh_triangle(xc, yc, 2*N, 2*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [0; 0; 0.2; 0.2];
yc = [0.2; 0.3; 0.4; 0.2];
dc = [0; 0.5; 0.25; 0];
[x1 y1 tri1 F1] = mesh_block_hermite(xc, yc, dc, 2*N, 2*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [0.2; 0.2; 1.2; 1.2];
yc = [0.2; 0.4; 0.3; 0.2];
dc = [0; 0.25; -0.25; 0];
[x1 y1 tri1 F1] = mesh_block_hermite(xc, yc, dc, 2*N, 10*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [1.2; 1.2; 1.6];
yc = [0.2; 0.3; 0.2];
[x1 y1 tri1 F1] = mesh_triangle(xc, yc, 2*N, 4*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);
% Presion
xc = [0; 0; -0.1];
yc = [0; -0.05; 0];
[x1 y1 tri1 F1] = mesh_triangle(xc, yc, N, 2*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [0; 0; 0.2; 0.2];
yc = [-0.05; 0; 0; -0.05];
[x1 y1 tri1 F1] = mesh_block(xc, yc, N, 2*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [0.2; 0.2; 1.2; 1.2];
yc = [-0.05; 0; 0; -0.05];
[x1 y1 tri1 F1] = mesh_block(xc, yc, N, 10*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

xc = [1.2; 1.2; 1.6];
yc = [0; -0.05; 0];
[x1 y1 tri1 F1] = mesh_triangle(xc, yc, N, 4*N);
[x y tri F] = mesh_meld(x, y, tri, F, x1, y1, tri1, F1, true);

dlmwrite('xx.txt', x, 'delimiter', ' ');
dlmwrite('yy.txt', y, 'delimiter', ' ');

% Término fuente
f = @(x,y) 0;

% Parámetros y condiciones de contorno
N = length(x); %Número de nodos
k = 20; % Conductividad del material compuesto [W/mK]
a = 0;  % Coeficeinte de u
b = 0;  % Coeficiente del grad
HTC = @(x,y) 1000;   % Coeficiente de transferencia de calor [W/m^2 K]
g = @(x,y) 293;      % [K] Temperatura del fluido interior
Text = @(x,y) 273;   % [K] Temperatura del fluido exterior

% Evaluar término independiente
F = zeros(N,1);
F(1:N) = f(x(1:N),y(1:N));

% Ensamblar matrices
[M R C D] = fem_mrc(x, y, tri, F, b, false);
A = a*M + b*C + k*R;

%Condición de contorno Dirichlet
front_D = load('frontd.txt');
for i=1:N
  for k=1:length(front_D)
    j = front_D(k);
    D(i) -= g(x(j), y(j))*A(i,j);
  end
end

A0 = A;
A0(front_D,:) = 0;
A0(:,front_D) = 0;

D(front_D) = g(x(front_D),y(front_D));
for j=1:length(front_D)
  i=front_D(j);
  A0(i,i) = 1;
end

%Condición de contorno Robin
front_R = load('frontr.txt');
G = ALPHA = sparse(N, 1);
ALPHA(front_R) = HTC(x(front_R),y(front_R));
G(front_R) = HTC(x(front_R),y(front_R)).*Text(x(front_R),y(front_R));
[AR DR] = fem_robin(x, y, tri, G, ALPHA);
D += DR;
A0 += AR;

% Resolución del sistema lineal del problema débil
uh = A0\D;

dlmwrite('uh5.txt', uh, 'delimiter', ' ');

%obtención de mínimas temperaturas
minT1=find(uh<=273.01);
minT2=find(uh<=273.005);
minimo=min(uh);
minT3=find(uh==minimo);


%Pintar mallado
figure (1)
triplot(tri, x, y,'-b');
hold on
plot(x(front_D),y(front_D),'-r','Linewidth',2)
plot(x(front_R),y(front_R),'-g','Linewidth',2)
title('Mallado junto con sus fronteras N=5')
legend('Mallado','Frontera interior','Frontera exterior')
grid minor

% Pintar Soución
figure (2)
trisurf(tri,x,y,uh);
shading("interp");
colormap("hot");
colorbar();
xlabel ("x");
ylabel ("y");
zlabel ("u");
title('Distribución de temperaturas del perfil alar N=5')
grid minor

%Pintar zona de mínima temperaturas
figure (3)
triplot(tri, x, y,'-b');
hold on
plot(x(minT1),y(minT1),'o','color','y','Linewidth',1.5)
plot(x(minT2),y(minT2),'o','color','g','Linewidth',1.5)
plot(x(minT3),y(minT3),'o','color','r','Linewidth',1.5)

title('Zona de mínima temperatura N=5')
legend('Mallado','Zona T<=273.01','Zona T<=273.005','Mínima temperatura')
grid minor

ind5=[109
190
1414
1412
1
10
702
1468
1470
91
100
827
1712
1714
191
200
1316
381
390
571
580
693];

figure(4)
triplot(tri, x, y,'-b');
hold on
plot(x(ind5),y(ind5),'o','color','r','Linewidth',2)

title('Puntos coincidentes N=5')
legend('Mallado','Puntos coincidentes')
grid minor

