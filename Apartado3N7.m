clear all;
clear;
clc;

% Depósito de combustible en perfil alar. Apartado 1

% Perfil %
% Fuente y borde de ataque
N = 7; % Densidad de mallado
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
% Presión
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

% Parámetros
N = length(x);         %Número de nodos
a = 0;                 %coeficiente de u
k = 3.8e-5;            %coeficiente de difusividad térmica
K = 20;                %conductividad térmica
b = 0;                 %coeficiente del grad u
F = zeros(N,1);        % Fuente
HTC = @(x,y) 1000;     %coeficiente de convección
Text = @(x,y) 273;     %temperatura exterior
neu = @(x,y) 0;

% Pasos temporales
dt = 5; % Delta de tiempo
T = 0:dt:250; %intervalo de tiempo
Nt = length(T); % Número de instantes temporales

% Evaluar término independiente
F = zeros(N,1);
F(1:N) = f(x(1:N),y(1:N));

%Condicion de contorno Neumann
front_N=load('frontn.txt');

% Filas -> nodos; Columnas -> tiempos

uht = zeros(N,Nt);

% Solución inicial
uht(:,1) = load('uh7.txt');

% Obtencion de las matrices
[M R C D] = fem_mrc(x, y, tri, F, b, false);
A = (1 + a*dt/2)*M + (dt/2)*(b*C+k*R);
B = (1 - a*dt/2)*M - (dt/2)*(b*C+k*R);

% Resolver
for t=2:Nt

  A0=A;

  printf("Resolviendo instante t=%f s\n", T(t));
  Bu = B*uht(:,t-1);

  %Condicion de contorno Neumann
  G = ALPHA = sparse(N, 1);
  G(front_N) = K * neu(x(front_N),y(front_N));
  [ARN DNN] = fem_robin(x, y, tri, G, ALPHA);
  Bu += DNN;

  %Condición de contorno Robin
  front_R = load('frontr.txt');
  G = ALPHA = sparse(N, 1);
  ALPHA(front_R) = HTC(x(front_R),y(front_R));
  G(front_R) = HTC(x(front_R),y(front_R)).*Text(x(front_R),y(front_R));
  [AR DR] = fem_robin(x, y, tri, G, ALPHA);
  A0 += AR;
  Bu += DR;
  uht(:,t) = A0\Bu;

endfor

% Pintar solución
figure(1);
colormap("hot");
for t=1:Nt
  trisurf(tri,x,y,uht(:,t));
  shading("interp");
  colorbar();
  caxis([273 293]); % fijar los límites de los ejes de color
  title (sprintf("Distribución de temperaturas al transcurso de t = %f s", T(t)));
  xlabel ("x");
  ylabel ("y");
  zlabel ("T");
  zlim([272 max(uht(:,t))]);
  view(0,90);
  grid minor
  %waitforbuttonpress();
   pause(0.1);
endfor

maxt=find(uht(:,t)==max(uht(:,t)));
max(uht(:,t))
maxT=maxt;

mintt=find(uht(:,t)<=273.005);
min(uht(:,t))
minT=mintt;

minTTd=find(uht(:,t)<=273.0001);
minTT=minTTd;

minTTT=find(uht(:,t)==273);

figure(2)
triplot(tri, x, y,'-g');
hold on
plot(x(maxT),y(maxT),'o','color','r','Linewidth',2)
plot(x(minT),y(minT),'o','color','b','Linewidth',2)
plot(x(minTT),y(minTT),'o','color','k','Linewidth',2)
plot(x(minTTT),y(minTTT),'o','color','c','Linewidth',2)
title(sprintf("Puntos térmicos clave del perfil para el tiemto: t = %f s", T(t)))
legend('Mallado','Punto de máxima temratura del perfil','Puntos con una T<=273.005','Puntos con una T<=273.0001','Punto de mínima temperatura')
grid minor

figure(3)
m=[1 50 125 250 350 750 1125 1500];
n=[292.60 288.94 286.48 283.45 281.09 276.70 274.54 273.65];
plot(m,n,'o')
title('Ajuste de curva a puntos de máxima temperatura del perfil del ala')
legend('Datos')
grid minor
